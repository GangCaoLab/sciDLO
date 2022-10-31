use std::collections::HashMap;
use std::fs::File;
use std::hash::Hash;
use std::io::{self, prelude::*};
use std::ops::AddAssign;
use std::path::Path;

extern crate flate2;
use flate2::read::GzDecoder;

pub fn open_file(path: &str) -> Box<dyn Read + Send + Sync> {
    let f: Box<dyn Read + Send + Sync> = if path.ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(path).unwrap()))
    } else {
        Box::new(File::open(path).unwrap())
    };
    f
}

pub fn open_file_buffered(path: &str) -> io::BufReader<Box<dyn Read + Send + Sync>> {
    let f = open_file(path);
    let buf_read = io::BufReader::new(f);
    buf_read
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn add_hashmap<T, R>(m1: HashMap<T, R>, m2: HashMap<T, R>) -> HashMap<T, R>
where
    T: Clone + Eq + Hash,
    R: Clone + AddAssign + Copy,
{
    let mut m = m1.clone();
    for (k, v) in m2 {
        *m.entry(k).or_insert(v) += v;
    }
    return m;
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
