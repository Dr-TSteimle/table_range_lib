use std::{
    collections::HashMap,
    fs::{File, rename},
    io::{BufRead, copy, ErrorKind, Write, BufReader},
    fmt::Debug,
    num::{ParseIntError, NonZeroUsize},
    path::Path
};

use flate2::read::GzDecoder;
use noodles::bgzf::{VirtualPosition, Reader as BGZReader, MultithreadedWriter};

#[derive(Debug)]
pub enum TRError {
    IoErr(std::io::Error),
    ParseIntError(ParseIntError),
    ColumnsNotFound,
    Infallible(std::convert::Infallible),
    InvalidFileFormat,
    NotSorted
}

impl From<std::io::Error> for TRError {
    fn from(e: std::io::Error) -> TRError {
        Self::IoErr(e)
    }
}

impl From<ParseIntError> for TRError {
    fn from(e: ParseIntError) -> Self {
        Self::ParseIntError(e)
    }
}

impl From<std::convert::Infallible> for TRError {
    fn from(e: std::convert::Infallible) -> Self {
        Self::Infallible(e)
    }
}

impl std::fmt::Display for TRError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TRError::IoErr(err) => write!(f, "[Error] Input/Output: {}", err.to_string()),
            TRError::ParseIntError(err) => write!(f, "ParseIntError {:?}", err),
            TRError::ColumnsNotFound => write!(f, "[Error] Error while parsing columns, are the arguments of the columns positions or the separator well filled ?"),
            TRError::Infallible(_err) => write!(f, "[Error] Error while parsing columns, are the arguments of the columns positions or the separator well filled ?"),
            TRError::InvalidFileFormat => write!(f, "[Error] Incompatible file format."),
            TRError::NotSorted => write!(f, "[Error] Positions aren't sorted."),
        }
    }
}

#[derive(Debug, Clone)]
struct GenomicPosition {
    contig: String,
    start : i32,
    end   : i32
}

impl GenomicPosition {
    fn from_line(line: &String, sep: &str, position_columns: &Vec<usize>) -> Result<GenomicPosition, TRError> {
        let str = line
        .split(sep)
        .enumerate()
        .filter(|(i, _)| position_columns.contains(i))
        .map(|(_,e)| e.to_owned())
        .collect::<Vec<String>>();

        if str.len() != 3 {
            return Err(TRError::ColumnsNotFound)
        }

        Ok(GenomicPosition {
            contig: str[0].clone(),
            start : str[1].trim().parse()?,
            end   : str[2].trim().parse()?,
        })
    }
}

struct  GenomicPositions(Vec<(u64, GenomicPosition, u64)>);

impl GenomicPositions {
    fn new_empty() -> GenomicPositions {
        GenomicPositions(Vec::new())
    }

    /// Load GenomicPositions from index file
    fn new_from_index(path: &str, sep: &str) -> Result<GenomicPositions, TRError> {
        let mut res = GenomicPositions(Vec::new());
        let mut reader = std::io::BufReader::new(File::open(&path)?);
        // let mut reader = File::open(&path).map(Reader::new)?;

        let mut n_line = 1;
        let mut line_buffer = String::new();
        loop {
            let line = reader.read_line(&mut line_buffer)?;
            if line == 0 {break;}
            
            let cols = line_buffer.split(sep).collect::<Vec<&str>>();
            if let (Ok(pos) , Ok(start), Ok(end), Ok(n_pos)) = (cols[0].trim().parse(), cols[2].trim().parse(), cols[3].trim().parse(), cols[4].trim().parse()) {
                res.push((
                    pos,
                    GenomicPosition { 
                        contig: cols[1].to_owned(),
                        start, end
                    },
                    n_pos
                ));
            } else {
                eprintln!("Error parsing index line: {} {:?}", n_line, cols);
                break;
            }

            line_buffer.clear();
            
            if n_line % 100_000 == 0 {
                println!("{} parsed", n_line);
            }
            n_line += 1;
        }
        println!("Index parsed {} ranges.", n_line);

        Ok(res)
    }

    fn push_from_line(&mut self, offset: u64, line: &String, sep: &str, position_columns: &Vec<usize>, n_pos: u64) -> Result<GenomicPosition, TRError> {
        let to_add = GenomicPosition::from_line(line, sep, position_columns)?;
        self.0.push((offset, to_add.clone(), n_pos));
        Ok(to_add)
    }

    fn write_index(&self, path: &str, by: usize) -> Result<(), TRError> {
        let mut write_buffer = std::io::BufWriter::new(File::create(&path)?);

        let mut last_contig = self.0[0].1.contig.clone();
        let mut pos_first = 0;
        
        for (i, current) in self.0.iter().enumerate() {
            if (i - pos_first) == by + 1 || last_contig != current.1.contig {
                let first = self.0.get(pos_first).unwrap();
                let last = self.0.get(i - 1).unwrap();
                let mut l = vec![
                    first.0.to_string(),
                    first.1.contig.clone(),
                    first.1.start.to_string(),
                    last.1.end.to_string(),
                    (i - 1 - pos_first).to_string()
                ].join("\t");
                l.push('\n');
                
                write_buffer.write(l.as_bytes())?;
                
                last_contig = current.1.contig.clone();
                pos_first = i;
            }
        }
        
        Ok(())
    }

    fn push(&mut self, value: (u64, GenomicPosition, u64)) {
        self.0.push(value)
    }

    fn positions(&self, positions: Vec<(String, i32)>, tolerance: i32) -> Vec<Option<(u64, u64)>> {
        let mut res = Vec::new();
        let mut positions_iter = positions.iter();

        let (current, mut next) = (positions_iter.next(), positions_iter.next());
        
        if let Some((contig, position)) = current {
            let mut current_contig = contig.to_string();
            let mut current_pos = *position;

            'A: for (offset, gp, n_pos) in self.0.iter() {
                if current_contig == gp.contig {
                    let mut to_next = false;
                    'B: loop {
                        if current_contig == gp.contig && gp.start - tolerance <= current_pos && current_pos <= gp.end + tolerance {
                            res.push(Some((*offset, *n_pos)));
                            to_next = true;
                        } else if let Some((next_contig, next_position)) = next {
                            if gp.contig == next_contig.to_string() {
                                if gp.start - tolerance <= *next_position && *next_position <= gp.end + tolerance {
                                    res.push(None);
                                    to_next = true;
                                }
                            }
                        }
    
                        if to_next {
                            if let Some((c, p)) = next {
                                (current_contig, current_pos) = (c.to_string(), *p);
                                next = positions_iter.next();
                            } else {
                                break 'A;
                            }
                        } else {
                            break 'B;
                        }
                    }
                }
                
            }
        }

        res
    }

    fn contigs_order(&self) -> Vec<String> {
        self.0.iter().map(|(_, gp, _)| gp.contig.to_string()).collect()
    }
}

pub struct TableFile {
    index: GenomicPositions,
    reader: BGZReader<File>
}

impl TableFile {
    pub fn new(path: &str, sep: &str, position_columns: &Vec<usize>, comment: &str) -> Result<TableFile, TRError> {
        let (index, reader) = TableFile::get_readers(path, sep, position_columns, comment)?;
        Ok(TableFile { index, reader })
    }

    fn get_readers(path: &str, sep: &str, position_columns: &Vec<usize>, comment: &str) -> Result<(GenomicPositions, BGZReader<File>), TRError> {
        let mut path = path.to_string();
        println!("Loading {}", path);

        if let Some(r) = Path::new(&path).extension() {
            if let Some(ext) = r.to_str() {
                if ext != "gz" {
                    let tp = format!("{}.gz", path.clone()).to_string();
                    if std::path::Path::new(&tp).exists() {
                        path = tp;
                    }
                }
            }
        }

        let mut reader = File::open(&path).map(BGZReader::new)?;
        let test = reader.seek(reader.virtual_position());
        let mut reader = match test {
            Ok(_) => reader,
            Err(err) => {
                match err.kind() {
                    ErrorKind::InvalidData => {
                        println!("The table file isn't in BGZ format.");
                        
                        let mut tmp_file = std::env::temp_dir();
                        tmp_file.push("reader-tmp.gz");
                        let mut writer = MultithreadedWriter::with_worker_count(NonZeroUsize::new(5).unwrap(), File::create(&tmp_file)?);
                        
                        match Path::new(&path).extension() {
                            Some(ext) => {
                                match ext.to_str() {
                                    Some("gz") => {
                                        println!("Converting from GZ format.");
                                        let mut decoder = GzDecoder::new(File::open(&path)?);
                                        copy(&mut decoder, &mut writer)?;
                                    },
                                    _ => {
                                        println!("Converting from TXT format.");
                                        let mut decoder = BufReader::new(File::open(&path)?);
                                        copy(&mut decoder, &mut writer)?;
                                        path = format!("{}.gz", path.clone()).to_string();
                                    },
                                }
                                
                            },
                            None => return Err(TRError::InvalidFileFormat),
                        }

                        rename(tmp_file, path.clone())?;
                        File::open(&path).map(BGZReader::new)?
                    },
                    _ => panic!("{:?}", err),
                }
            },
        };

        let id = format!("{}.id", path);
        let index_path = Path::new(&id);

        let mut gp = GenomicPositions::new_empty();
        if index_path.exists() {
            println!("Loading index...");
            gp = GenomicPositions::new_from_index(&id, "\t")?;
        } else {
            println!("Saving index...");

            let mut n_line = 0;
            let mut line_buffer = String::new();

            let mut last_gp: Option<GenomicPosition> = None;
            let n_chr_comment = comment.len();

            loop {
                let offset = reader.virtual_position().try_into().unwrap();

                match reader.read_line(&mut line_buffer) {
                    Ok(size) => {
                        if size == 0 || line_buffer.len() < n_chr_comment { break; }
                        if &line_buffer[..n_chr_comment] != comment { 
                            match gp.push_from_line(offset, &line_buffer, sep, &position_columns, 1) {
                                Ok(gp) => {
                                    if let Some(ls) = last_gp {
                                        if ls.start > gp.start && ls.contig == gp.contig {
                                            eprintln!("last start {}:{}, current start {}:{}", ls.contig, ls.start, gp.contig, gp.start);
                                            return Err(TRError::NotSorted);
                                        }
                                    }
                                    last_gp = Some(gp);
                                },
                                Err(_) => {
                                    println!("Parsing error line: {}\n{}", n_line, &line_buffer);
                                    return Err(TRError::ColumnsNotFound)
                                },
                            }
                        }
                    },
                    Err(err) => {
                        panic!("{:?}", err);
                    }
                }
                line_buffer.clear();

                n_line += 1;

                if n_line % 100_000 == 0 {
                    println!("{} parsed", n_line);
                }
            }

            gp.write_index(&id, 1_000)?;    
        }
        Ok((gp, reader))
    }

    pub fn get_positions_lines(&mut self, positions: Vec<(String, i32)>, sep: &str, position_columns: &Vec<usize>, tolerance: i32, comment: &str) -> Result<Vec<Option<String>>, TRError> {
        let ordered = order_positions(positions, self.index.contigs_order());
        // println!("ordered {:?}", ordered);
        let positions = self.index.positions(ordered.iter().map(|(_, p)| p.clone()).collect(), tolerance);
        // println!("positions {:?}", positions);

        let mut dedup: HashMap<u64, (Vec<usize>, u64)> = HashMap::new();
        for (i, p) in positions.iter().enumerate() {
            if let Some((pos, max_lines)) = p {
                if let Some(x) = dedup.get_mut(pos) {
                    x.0.push(i);
                } else {
                    dedup.insert(*pos, (vec![i], *max_lines));
                }
            }
        }
        let mut dedup = Vec::from_iter(dedup.iter());
        dedup.sort_by(|a,b| a.0.cmp(b.0));
        // println!("dedup {:?}", dedup);

        let mut res: Vec<(usize, Option<String>)> = ordered.iter().map(|(index, _)| (*index, None)).collect();
        let n_chr_comment = comment.len();

        for (offset, (items_ids, max_lines)) in dedup.into_iter() {
            let mut line_buffer = String::new();
            let mut n_lines = 0;

            let mut curr_i = 0;

            let mut ordered_id = *items_ids.get(curr_i).unwrap();
            let (_, (_, mut position)) = ordered.get(ordered_id).unwrap().clone();

            self.reader.seek(VirtualPosition::from(*offset)).unwrap();

            loop {
                match self.reader.read_line(&mut line_buffer) {
                    Ok(code) => {
                        if code == 0 { break; } else {
                            if n_lines == max_lines + 1 || line_buffer.len() < n_chr_comment { break; }

                            if &line_buffer[..n_chr_comment] != comment {

                                let str = line_buffer
                                .split(sep)
                                .enumerate()
                                .filter(|(i, _)| position_columns.contains(i))
                                .map(|(_,e)| e.to_owned())
                                .collect::<Vec<String>>();

                                if str.len() != 3 {
                                    return Err(TRError::ColumnsNotFound)
                                }
                                let (start, stop) = (str[1].trim().parse::<i32>()?, str[2].trim().parse::<i32>()?);

                                loop {
                                    if start - tolerance <= position && position <= stop + tolerance {
                                        let mut res_item = res.get_mut(ordered_id).unwrap();
                                        res_item.1 = Some(line_buffer.clone().trim().to_string()); 
        
                                        curr_i += 1;
                                        if let Some(o_i) = items_ids.get(curr_i) {
                                            ordered_id = *o_i;
                                            if ordered.get(ordered_id).is_some() {
                                                (_, (_, position)) = ordered.get(ordered_id).unwrap().clone();
                                            } else { break; }
                                        } else { break; }
                                    } else { break; }
                                }
                            }
                            
                            line_buffer.clear();
                        }
                    },
                    Err(_) => panic!("Error parsing file"),
                }
                n_lines += 1;
            }
        }

        res.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(res.iter().map(|(_,a)| a.clone()).collect())
    }
}

fn order_positions (input: Vec<(String, i32)>, contigs_order: Vec<String>) -> Vec<(usize, (String, i32))> {
    let mut input: Vec<(usize, (String, i32))> = input.into_iter().enumerate().collect();
    input.sort_by(|(_, (a_contig, a_pos)), (_, (b_contig, b_pos))| {
        if a_contig != b_contig {
            for c in contigs_order.iter() {
                if c == a_contig {
                    return std::cmp::Ordering::Less
                } else if c == b_contig {
                    return std::cmp::Ordering::Greater
                }
            }
            std::cmp::Ordering::Greater
        } else {
            a_pos.cmp(b_pos)
        }
    });

    input
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Instant;

    #[test]
    fn it_works() {
        let now = Instant::now();
    
        let path = "/home/thomas/NGS/ref/hg19/rmsk.txt.gz";
        let sep = "\t";
        let position_columns = vec![5, 6, 7];
        let comment = "#";

        let mut res = TableFile::new(path, sep, &position_columns, comment).unwrap();
        let my_pos = vec![("chr14".to_string(), 22555233),("chr1".to_string(), 249_240_620), ("chr10".to_string(), 524_779_845), ("chr1".to_string(), 249_240_621)];

        let expected = vec![
            None,
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string()), 
            None,
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string())
        ];

        if let Ok(res) = res.get_positions_lines(my_pos.clone(), sep, &position_columns, 0, comment) {
            for r in my_pos.into_iter().zip(&res) {
                println!("{:?}", r.1);
            }

            for r in res.into_iter().zip(expected) {
                assert_eq!(r.0, r.1);
            }
        }

        println!("{:#?}", now.elapsed());
        
    }

    #[test]
    fn it_works2() {
        let now = Instant::now();
    
        let path = "/home/thomas/NGS/ref/hg19/gencode.v28lift37.basic.annotation.gtf.gz";
        let sep = "\t";
        let position_columns = vec![0, 3, 4];
        let comment = "#";

        let mut res = TableFile::new(path, sep, &position_columns, comment).unwrap();
        let my_pos = vec![("chr14".to_string(), 19_013_295), ("chr14".to_string(), 105259757)];
        
        // let my_pos = vec![("chr14".to_string(), 105259757)];

        let expected = vec![
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string()), 
            None,
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string())
        ];

        if let Ok(res) = res.get_positions_lines(my_pos.clone(), sep, &position_columns, 0, comment) {
            for r in my_pos.into_iter().zip(&res) {
                println!("{:?}", r.1);
            }

            // for r in res.into_iter().zip(expected) {
            //     assert_eq!(r.0, r.1);
            // }
        }

        println!("{:#?}", now.elapsed());
        
    }
    
}
