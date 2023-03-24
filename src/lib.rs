use std::collections::HashMap;
use std::fs::{File, rename};
use std::io::{BufRead, copy, ErrorKind, Write};
use std::fmt::{Formatter, Debug, self};
use std::num::{ParseIntError, NonZeroUsize};
use std::path::Path;

use flate2::read::GzDecoder;
use noodles::bgzf::{VirtualPosition, Reader, MultithreadedWriter};

pub enum TRError {
    IoErr(std::io::Error),
    ParseIntError(ParseIntError),
    ColumnsNotFound,
    Infallible(std::convert::Infallible)
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

impl Debug for TRError {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        match self {
            TRError::IoErr(err) => write!(f, "[ERROR] Input/Output: {}", err.to_string()),
            TRError::ParseIntError(err) => write!(f, "ParseIntError {:?}", err),
            TRError::ColumnsNotFound => write!(f, "ColumnsNotFound"),
            TRError::Infallible(_err) => write!(f, "ColumnsNotFound"),
        }
    }
}

#[derive(Debug)]
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

    fn push_from_line(&mut self, offset: u64, line: &String, sep: &str, position_columns: &Vec<usize>, n_pos: u64) -> Result<(), TRError> {
        self.0.push((offset, GenomicPosition::from_line(line, sep, position_columns)?, n_pos));
        Ok(())
    }

    fn write_index(&self, path: &str, by: usize) -> Result<(), TRError> {
        let mut write_buffer = std::io::BufWriter::new(File::create(&path)?);

        let mut last_contig = self.0[0].1.contig.clone();
        let mut pos_first = 0;
        
        for (i, current) in self.0.iter().enumerate() {

            if (i - pos_first) == by + 1|| last_contig != current.1.contig {
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
        let mut curr_to_find = 0;
        let len_positions = positions.len();
        let mut has_search_in_contig = false;
        for (offset, gp, n_pos) in self.0.iter() {
            if curr_to_find == len_positions { break; }
            if gp.contig == positions[curr_to_find].0 {
                has_search_in_contig = true;

                if positions[curr_to_find].1 < gp.start { continue;}
                if positions[curr_to_find].1 > gp.end { continue;}

                loop {
                    if curr_to_find == len_positions { break; }
                    if gp.contig == positions[curr_to_find].0 && gp.start - tolerance <= positions[curr_to_find].1 && positions[curr_to_find].1 <= gp.end + tolerance {
                        res.push(Some((*offset, *n_pos)));
                        curr_to_find += 1;
                    } else { break; }
                }
            } else if has_search_in_contig {
                // after not found
                res.push(None);
                curr_to_find += 1;
                has_search_in_contig = false;
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
    reader: Reader<File>
}

impl TableFile {
    pub fn new(path: &str, sep: &str, position_columns: &Vec<usize>) -> Result<TableFile, TRError> {
        let (index, reader) = TableFile::get_readers(path, sep, position_columns)?;
        Ok(TableFile { index, reader })
    }

    fn get_readers(path: &str, sep: &str, position_columns: &Vec<usize>) -> Result<(GenomicPositions, Reader<File>), TRError> {
        let mut reader = File::open(&path).map(Reader::new)?;
        
        let test = reader.seek(reader.virtual_position());

        let mut reader = match test {
            Ok(_) => reader,
            Err(err) => {
                match err.kind() {
                    ErrorKind::InvalidData => {
                        println!("Not BGZF, converting...");
                        
                        let mut tmp_file = std::env::temp_dir();
                        tmp_file.push("reader-tmp.gz");
                        
                        let mut decoder = GzDecoder::new(File::open(path)?);

                        let mut writer = MultithreadedWriter::with_worker_count(NonZeroUsize::new(5).unwrap(), File::create(&tmp_file)?);
                        copy(&mut decoder, &mut writer)?;

                        rename(tmp_file, path)?;
                        File::open(&path).map(Reader::new)?
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

            loop {
                let offset = reader.virtual_position().try_into().unwrap();

                match reader.read_line(&mut line_buffer) {
                    Ok(size) => {
                        if size == 0 { break; }
                        
                        match gp.push_from_line(offset, &line_buffer, sep, &position_columns, 1) {
                            Ok(_) => (),
                            Err(_) => println!("Parsing error line: {}\n{}", n_line, &line_buffer),
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

    pub fn get_positions_lines(&mut self, positions: Vec<(String, i32)>, sep: &str, position_columns: &Vec<usize>, tolerance: i32) -> Result<Vec<Option<String>>, TRError> {
        let ordered = order_positions(positions, self.index.contigs_order());
        let positions = self.index.positions(ordered.iter().map(|(_, p)| p.clone()).collect(), tolerance);
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

        let mut res: Vec<(usize, Option<String>)> = ordered.iter().map(|(index, _)| (*index, None)).collect();

        for (offset, (items_ids, max_lines)) in dedup.into_iter() {
            let mut line_buffer = String::new();
            let mut n_lines = 0;

            let mut curr_i = 0;

            // let (mut index, (mut contig, mut position)) = sub_ordered.get(curr_i).unwrap().clone();
            let mut ordered_id = *items_ids.get(curr_i).unwrap();
            let (_, (_, mut position)) = ordered.get(ordered_id).unwrap().clone();

            self.reader.seek(VirtualPosition::from(*offset)).unwrap();

            loop {
                match self.reader.read_line(&mut line_buffer) {
                    Ok(code) => {
                        if code == 0 { break; } else {
                            if n_lines == max_lines + 1 { break; }

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
                                        (_, (_, position)) = ordered.get(ordered_id).unwrap().clone();
                                    } else {
                                        break;
                                    }
                                } else {
                                    break;
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

        let mut res = TableFile::new(path, sep, &position_columns).unwrap();
        let my_pos = vec![("chr1".to_string(), 249_240_620), ("chr10".to_string(), 524_779_845), ("chr1".to_string(), 249_240_621)];

        let expected = vec![
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string()), 
            None,
            Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string())
        ];

        if let Ok(res) = res.get_positions_lines(my_pos.clone(), sep, &position_columns, 0) {
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
    
        let path = "/home/thomas/NGS/ref/hg19.refGene.sortedhh.gtf.gz";

        let sep = "\t";
        let position_columns = vec![0, 3, 4];

        let mut res = TableFile::new(path, sep, &position_columns).unwrap();
        let my_pos = vec![("chr1".to_string(), 249_240_620), ("chr10".to_string(), 524_779_845), ("chr1".to_string(), 249_240_621)];

        // let expected = vec![
        //     Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string()), 
        //     None,
        //     Some("2486\t1442\t13\t8\t38\tchr1\t249239883\t249240621\t-10000\t+\t(TTAGGG)n\tSimple_repeat\tSimple_repeat\t1\t716\t0\t3".to_string())
        // ];

        if let Ok(res) = res.get_positions_lines(my_pos.clone(), sep, &position_columns, 0) {
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
