use std::{io::{self, BufRead, ErrorKind, copy, BufReader}, path::Path, fs::{File, rename}, num::NonZeroUsize, ops::RangeInclusive, collections::HashMap};

use flate2::read::GzDecoder;
use noodles_core::{Position, region::Interval};
use noodles_bgzf as bgzf;
use bgzf::Reader;
use noodles_csi::{index::reference_sequence::bin::Chunk, BinningIndex};
use noodles_tabix as tabix;
use tabix::Index;
use rangemap::RangeInclusiveMap;
use rayon::prelude::*;

fn parse_record(s: &str, separator: String, columns_positions: Vec<usize>) -> io::Result<(&str, i32, i32)> {
    let components = s
    .split(&separator)
    .collect::<Vec<&str>>();
        
    let components = components.into_iter()
        .enumerate()
        .filter(|(i, _)| columns_positions.contains(i))
        .map(|(_, e)| e)
        .collect::<Vec<&str>>();

    let mut components = components.iter();

    let reference_sequence_name = components
        .next()
        .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?;

    let start = components
        .next()
        .and_then(|t|if t == &"0" {"1".parse().ok()} else {t.parse().ok()})
        .ok_or_else(|| {
            println!("start {:?}", s);
            io::Error::from(io::ErrorKind::InvalidData)
        })?;

    let end = components
        .next()
        .and_then(|t| t.parse().ok())
        .ok_or_else(|| {
            println!("start {:?}", s);
            io::Error::from(io::ErrorKind::InvalidData)
        })?;

    Ok((reference_sequence_name, start, end))
}

pub fn get_readers(path: &str, seperator: String, columns_positions: Vec<usize>, comment: &str) -> io::Result<(tabix::Index, bgzf::Reader<File>)> {
    let mut path = path.to_string();
    println!("Loading {}", path);

    // if file extension is not gz but a gz file exists add .gz to the path
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
    
    let mut reader = File::open(&path).map(bgzf::Reader::new).unwrap();
    
    // If gzip and not bgzip convert using flate2 reader and noodles bgzip writer
    let test = reader.seek(reader.virtual_position());
    let mut reader = match test {
        Ok(_) => reader,
        Err(err) => {
            match err.kind() {
                ErrorKind::InvalidData => {
                    println!("The table file isn't in BGZ format.");
                    
                    let mut tmp_file = std::env::temp_dir();
                    tmp_file.push("reader-tmp.gz");

                    let mut writer = bgzf::MultithreadedWriter::with_worker_count(NonZeroUsize::new(2).unwrap(), File::create(&tmp_file)?);
                    
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
                        None => return Err(io::Error::from(io::ErrorKind::InvalidData)),
                    }

                    rename(tmp_file, path.clone())?;
                    File::open(&path).map(bgzf::Reader::new)?
                },
                _ => panic!("{:?}", err),
            }
        },
    };

    // load the index
    let id = format!("{}.tbi", path);
    let index_path = Path::new(&id);

    if !index_path.exists() {
        println!("Saving index...");

        let n_chr_comment = comment.len();

        let mut indexer = tabix::Index::indexer();
        indexer.set_header(tabix::index::header::Builder::bed().build());

        let mut buf = String::new();
        let mut start_position = reader.virtual_position();

        loop {
            buf.clear();

            if reader.read_line(&mut buf)? == 0 {
                break;
            }

            let end_position = reader.virtual_position();
            let chunk = Chunk::new(start_position, end_position);
            
            if buf.len() <= n_chr_comment { continue; } else {
                if &buf[..n_chr_comment] == comment { continue; }
            }

            let (reference_sequence_name, start, end) = parse_record(buf.trim_end(), seperator.to_string(), columns_positions.to_vec())?;
            indexer.add_record(reference_sequence_name, Position::new(start as usize).unwrap(), Position::new(end as usize).unwrap(), chunk);

            start_position = end_position;
        }

        let index = indexer.build();
        let mut writer = tabix::Writer::new(File::create(&index_path)?);
        writer.write_index(&index)?;
    }

    let index = tabix::read(index_path)?;
    Ok((index, reader))
}

pub fn get_position(start: usize, index: Index, reader: &mut bgzf::Reader<File>, separator: String, columns_positions: Vec<usize>) -> RangeInclusiveMap::<i32, String> {
    let pos = Interval::from(RangeInclusive::new(
        Position::new(start).unwrap(),
        Position::new(start + 1).unwrap()
    ));

    let query = index.query(0, pos).unwrap();
    let mut rangemap = RangeInclusiveMap::new();
    for chunk in query.iter() {
        let _ = reader.seek(chunk.start());
        let mut buf = String::new();
        loop {
            buf.clear();

            if reader.read_line(&mut buf).unwrap() == 0 { break; }

            let (_reference_name, start, end) = parse_record(&buf, separator.to_string(), columns_positions.to_vec()).unwrap();
            rangemap.insert(start..=end, buf.clone());

            if reader.virtual_position() == chunk.end() { break; }
        }
    }
    rangemap
}

pub struct TableFile {
    pub reader: Reader<File>,
    pub index: Index,
    pub options: TableFileOpts
}

#[derive(Clone)]
pub struct TableFileOpts {
    pub path: String,
    pub separator: String,
    pub comment  : String,
    pub position_columns: Vec<usize>,
    pub tolerance: i32,
}

impl TableFileOpts {
    pub fn new(path: String, separator: String, comment: &str, position_columns: Vec<usize>, tolerance: i32) -> TableFileOpts {
        TableFileOpts {
            path: path.to_string(),
            separator,
            comment: comment.to_string(),
            position_columns,
            tolerance
        }
    }
}

impl TableFile {
    pub fn open(options: TableFileOpts) -> io::Result<TableFile> {
        let (index, reader) = get_readers(&options.path, options.separator.to_string(), options.position_columns.to_vec(), &options.comment)?;
        Ok(TableFile { reader, index, options })
    }

    pub fn get_positions(&mut self, positions: Vec<(String, i32)>) -> io::Result<Vec<Option<Vec<String>>>> {
        let pos_len = positions.len();
        let positions_id: Vec<(usize, Option<(usize, i32)>)> = positions.into_iter().enumerate().map(|(i, (c, p) )| {
            if let Some(ci) = self.get_reference_id(&c) { (i, Some((ci, p))) } else { (i, None) }}
        ).collect();

        // group positions (values) by reference name (key) in an Hashmap
        let mut by_ref: HashMap<usize, Vec<(usize, i32)>> = HashMap::new();
        for (id, cp) in positions_id.iter() {
            if let Some((c, p)) = cp {
                if let Some(c_positions) = by_ref.get_mut(c) {
                    c_positions.push((*id, *p));
                } else {
                    by_ref.insert(*c, vec![(*id, *p)]);
                }
            }
        }

        let mut positions_acc: Vec<Option<Vec<String>>> = vec![None; pos_len];
        // for each ref
        for (ci, mut positions) in by_ref.into_iter() {
            positions.sort_by(|a, b| a.1.cmp(&b.1));
            
            let mut chunks: Vec<Chunk> = Vec::new();
            // for each position of the ref retrieve corresponding chunks
            for (_, position) in positions.iter() {
                let position = Position::try_from(*position as usize).unwrap();
                let interval = position..=position;
                chunks.extend(self.index.query(ci, interval).unwrap());
            }

            // read and parse lines of chunk and add to RangeInclusiveMap
            let mut inc_map = RangeInclusiveMap::new();
            for chunk in chunks.iter() {
                let _ = self.reader.seek(chunk.start());
                let mut buf = String::new();
                loop {
                    buf.clear();
        
                    if self.reader.read_line(&mut buf).unwrap() == 0 { break; }
        
                    let (_reference_name, start, end) = parse_record(&buf, self.options.separator.to_string(), self.options.position_columns.to_vec()).unwrap();
                    inc_map.insert((start - self.options.tolerance)..=(end +  self.options.tolerance), buf.clone());
        
                    if self.reader.virtual_position() == chunk.end() { break; }
                }
            }

            // for each positions retieve overlapping intervals contained in RangeInclusiveMap
            for (id, p) in positions.iter() {
                let interval = *p..=*p;
                let r: Vec<String> = inc_map
                    .overlapping(&interval)
                    .map(|(_, s)| s.trim().to_string())
                    .collect();
                
                let e = positions_acc.get_mut(*id).unwrap();
                if r.len() > 0 {
                    *e = Some(r);
                } else {
                    *e = None;
                }
            }
        }

        Ok(positions_acc)
    }

    pub fn get_reference_id(&self, reference_name: &str) -> Option<usize> {
        self.index.header().reference_sequence_names().iter().enumerate().filter(|(_, e)| *e == reference_name).map(|(i, _)| i).nth(0)
    }
}

pub fn positions_par(positions: Vec<(String, i32)>, options: TableFileOpts, chunk_size: usize) -> io::Result<Vec<Option<Vec<String>>>> {
    let positions: Vec<(usize, (String, i32))> = positions.into_iter().enumerate().collect();
    type ItemVec = (usize, Option<Vec<String>>);
    let mut res: Vec<(usize, Option<Vec<String>>)> = positions
        .par_chunks(chunk_size)
        .map(|positions| -> io::Result<Vec<ItemVec>> {
            let mut table_file = TableFile::open(options.clone())?;
            let res = table_file.get_positions(positions.into_iter().map(|(_, (c, p))| (c.to_string(), *p)).collect())?;
            let res: Vec<ItemVec> = res.iter().zip(positions.iter()).map(|(r,(id, _))| (*id, r.clone())).collect();
            Ok(res)
        })
        .filter_map(|e| e.ok())
        .flat_map(|e| { e })
        .collect();
    res.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(res.into_iter().map(|(_, e)| e).collect())
}

#[cfg(test)]
mod tests {
    use std::time::Instant;
    use super::*;

    #[test]
    fn it_works() {
        let options = TableFileOpts::new(
            "data/rmsk.txt.gz".to_string(),
            '\t'.to_string(),
            "#",
            vec![5, 6, 7],
            0
        );
        
        let pos = vec![
            ("chr1".to_string(), 249_239_882),
            ("chr1".to_string(), 47_692_481),
            ("chr14".to_string(), 22_555_233),
            ("chr1".to_string(), 249_240_620),
            ("chr1".to_string(), 249_239_883),
            ("chr10".to_string(), 524_779_845),
            ("chr1".to_string(), 249_240_621),
            ("chr1".to_string(), 23_636_654),
        ];
        let all_pos: Vec<(String, i32)> = pos.into_iter().cycle().take(1000).collect();
        
        let now = Instant::now();
        let mut table_file = TableFile::open(options.clone()).unwrap();
        let _res = table_file.get_positions(all_pos.clone()).unwrap();
        println!("{}ms", now.elapsed().as_millis());

        let now = Instant::now();
        let _res = positions_par(all_pos, options, 200);
        println!("par {}ms", now.elapsed().as_millis());
        
        // res.iter().for_each(|e| println!("{:?}", e));
    }

    #[test]
    fn it_works2() {
        let options = TableFileOpts::new("/home/thomas/NGS/ref/hg19/gencode.v28lift37.basic.annotation.gtf.gz".to_string(), '\t'.to_string(), "#", vec![0, 3, 4], 0);
        
        let pos = vec![
            ("chr14".to_string(), 19_013_295),
            ("chr14".to_string(), 19_038_381),
            ("chr14".to_string(), 19_494_423),
            ("chr14".to_string(), 19_488_487),
            ("chr14".to_string(), 19_529_869),
            ("chr14".to_string(), 19_529_941),
            ("chr14".to_string(), 19_817_189),
        ];
        let all_pos: Vec<(String, i32)> = pos.into_iter().cycle().take(1000).collect();
        
        let now = Instant::now();
        let mut table_file = TableFile::open(options.clone()).unwrap();
        let res = table_file.get_positions(all_pos.clone()).unwrap();
        println!("{}ms", now.elapsed().as_millis());

        let now = Instant::now();
        let res = positions_par(all_pos.clone(), options.clone(), 200).unwrap();
        println!("par {}ms, n = {}", now.elapsed().as_millis(), res.len());

        // res.iter().for_each(|e| println!("-> {:?}", e));
    }
}
