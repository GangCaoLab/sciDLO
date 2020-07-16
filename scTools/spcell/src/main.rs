use std::fs::File;
use std::thread;
use std::sync::mpsc;
use std::sync::{Mutex, Arc};
use std::time::Duration;
use std::collections::HashMap;
use std::fmt;
use std::io::{self, prelude::*};
use std::ops;

extern crate bio;
extern crate clap;
extern crate log;
extern crate simple_logger;
extern crate strsim;
extern crate vpsearch;

use clap::{Arg, App};
use log::{info};
use strsim::{hamming};
use regex::Regex;
use bio::alphabets::dna::revcomp;

use utils::{read_lines, open_file, add_hashmap};



fn load_barcodes(barcodes_file: &str) -> Vec<String> {
    let mut codes = vec![];
    if let Ok(lines) = read_lines(barcodes_file) {
        for line in lines {
            let line = line.unwrap();
            let code = line.trim_end().to_string();
            codes.push(code);
        }
    }
    codes
}


struct BarcodeSearch {
    codes: Vec<String>,
    dist_thresh: usize,
}

impl BarcodeSearch {
    fn new(codes: Vec<String>, dist_thresh: usize) -> Self {
        Self {
            codes: codes,
            dist_thresh: dist_thresh,
        }
    }

    fn search(&self, code: &String) -> Option<(&String, usize)> {
        let mut min_dist = code.len();
        let mut index: usize = 0;
        for (i, c) in self.codes.iter().enumerate() {
            let dist = hamming(c, code).unwrap();
            if dist < min_dist {
                min_dist = dist;
                index = i;
            }
        }
        let res = if min_dist > self.dist_thresh {
            None
        } else {
            Some((&self.codes[index], min_dist))
        };
        res
    }
}

struct PairRec {
    code_r1_l: String,
    code_r2_l: String,
    code_r1_r: String,
    code_r2_r: String,
}


impl PairRec {
    fn from_line(line: &str, re_codes: &Regex) -> Self {
        let line = line.trim_end();
        let codes = re_codes.captures(line).unwrap();
        let r1_r_rc = String::from_utf8(revcomp(codes[3].as_bytes().to_vec())).unwrap();
        let r2_r_rc = String::from_utf8(revcomp(codes[4].as_bytes().to_vec())).unwrap();
        Self {
            code_r1_l: codes[1].to_string(),
            code_r2_l: codes[2].to_string(),
            code_r1_r: r1_r_rc,
            code_r2_r: r2_r_rc,
        }
    }
}

#[derive(Clone)]
struct Counter {
    b1b2_not_match: u64,
    r1r2_not_match: u64,
    barcode_not_found: u64,
    valid: u64,
    total: u64,
    barcode_cnts: HashMap<String, u64>,
}

impl Counter {
    fn new() -> Self {
        Self {
            b1b2_not_match: 0,
            r1r2_not_match: 0,
            barcode_not_found: 0,
            valid: 0,
            total: 0,
            barcode_cnts: HashMap::new(),
        }
    }
}

impl ops::Add<Counter> for Counter {
     type Output = Counter;

     fn add(self, _rhs: Counter) -> Counter {
          Self {
            b1b2_not_match: self.b1b2_not_match + _rhs.b1b2_not_match,
            r1r2_not_match: self.r1r2_not_match + _rhs.r1r2_not_match,
            barcode_not_found: self.barcode_not_found + _rhs.barcode_not_found,
            valid: self.valid + _rhs.valid,
            total: self.total + _rhs.total,
            barcode_cnts: add_hashmap(self.barcode_cnts, _rhs.barcode_cnts),
          }
     }
}

impl fmt::Display for Counter {
     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
          let ratio = |c| {
               if self.total == 0 { return format!("0%"); }
               format!("{:.2}%", ((c*100) as f64) / (self.total as f64))
          };
          let mut msg = format!("Count result:\n\
               valid:\t{}\t{}\n\
               b1b2_not_match\t{}\t{}\n\
               r1r2_not_match\t{}\t{}\n\
               barcode_not_found\t{}\t{}\n\
               total reads\t{}\n\n",
               self.valid, ratio(self.valid),
               self.b1b2_not_match, ratio(self.b1b2_not_match),
               self.r1r2_not_match, ratio(self.r1r2_not_match),
               self.barcode_not_found, ratio(self.barcode_not_found),
               self.total,
          );
          msg.push_str("Barcode counts:\n");
          let mut items_bar_cnts: Vec<(&String, &u64)> = self.barcode_cnts.iter().collect();
          items_bar_cnts.sort_by(|t1, t2| t2.1.cmp(t1.1));
          for (barcode, cnt) in items_bar_cnts {
               msg.push_str(&format!("{}\t{}\n", barcode, cnt));
          }
          write!(f, "{}", msg)
     }
}



fn locate_barcode(
        pair_rec: &PairRec,
        search: &BarcodeSearch,
        max_diff_b1b2: usize,
        max_diff_r1r2: usize,
        counter: &mut Counter,
        ) -> Option<String> {
    counter.total += 1;
    if hamming(&pair_rec.code_r1_l, &pair_rec.code_r1_r).unwrap() > max_diff_b1b2 ||
       hamming(&pair_rec.code_r2_l, &pair_rec.code_r2_r).unwrap() > max_diff_b1b2
    {
        counter.b1b2_not_match += 1;
        return None;
    }
    if hamming(&pair_rec.code_r1_l, &pair_rec.code_r2_l).unwrap() > max_diff_r1r2 ||
       hamming(&pair_rec.code_r1_r, &pair_rec.code_r2_r).unwrap() > max_diff_r1r2
    {
        counter.r1r2_not_match += 1;
        return None 
    }

    let codes = [&pair_rec.code_r1_l, &pair_rec.code_r1_r,
                 &pair_rec.code_r2_l, &pair_rec.code_r2_r];
    let mut res_vec = Vec::with_capacity(4);
    for code in &codes {
        if let Some((barcode, dist)) = search.search(code) {
            res_vec.push((barcode, dist))
        }
    }
    if res_vec.len() > 0 {
        counter.valid += 1;
        res_vec.sort_by(|t1, t2| t1.1.partial_cmp(&t2.1).unwrap());
        let code = res_vec[0].0.clone();
        *counter.barcode_cnts.entry(code.clone()).or_insert(0) += 1;
        return Some(code)
    } else {
        counter.barcode_not_found += 1;
        return None
    }
}


fn main() {
    simple_logger::init().unwrap();

    let matches = App::new("Split pairs file by barcodes.")
        .arg(Arg::with_name("pairs_file")
            .required(true)
            .help("Path of input pairs file."))
        .arg(Arg::with_name("barcodes_file")
            .required(true)
            .help("Path to the txt file which store all barcodes."))
        .arg(Arg::with_name("dist_thresh")
            .short("d")
            .long("dist_thresh")
            .takes_value(true)
            .help("Threshould of edit distance between reads \
                   barcodes to the barcodes in provided file."))
        .arg(Arg::with_name("max_diff_b1b2")
            .short("m")
            .long("max_diff_b1b2")
            .takes_value(true)
            .help("Threshould of edit distance between barcodes in same linker."))
        .arg(Arg::with_name("max_diff_r1r2")
            .short("M")
            .long("max_diff_r1r2")
            .takes_value(true)
            .help("Threshould of edit distance between R1 barcode to R2 barcode."))
        .arg(Arg::with_name("output_prefix")
            .short("o")
            .long("output_prefix")
            .required(true)
            .takes_value(true)
            .help("Prefix of output files."))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("Number of threads used for processing reads."))
        .arg(Arg::with_name("wait_timeout")
            .long("wait_timeout")
            .takes_value(true)
            .help("Wait time for end channel timeout."))
        .get_matches();

    let pairs_path = matches.value_of("pairs_file").unwrap();
    let barcodes_path = matches.value_of("barcodes_file").unwrap();
    let dist_thresh = matches.value_of("dist_thresh").unwrap_or("1");
    let dist_thresh: usize = dist_thresh.parse().unwrap();
    let max_diff_b1b2 = matches.value_of("max_diff_b1b2").unwrap_or("2");
    let max_diff_b1b2: usize = max_diff_b1b2.parse().unwrap();
    let max_diff_r1r2 = matches.value_of("max_diff_r1r2").unwrap_or("4");
    let max_diff_r1r2: usize = max_diff_r1r2.parse().unwrap();
    let output_prefix = matches.value_of("output_prefix").unwrap();
    let threads = matches.value_of("threads").unwrap_or("1");
    let threads: u8 = threads.parse().unwrap();
    let wait_t = matches.value_of("wait_timeout").unwrap_or("500");
    let wait_t: u64 = wait_t.parse().unwrap();

    info!("pairs_file: {} barcodes_file: {}\n \
           dist_thresh: {}, max_diff_b1b2: {}, max_diff_r1r2: {}\n \
           threads: {} wait_timeout: {}",
           pairs_path, barcodes_path, dist_thresh, max_diff_b1b2, max_diff_r1r2, threads, wait_t);

    let barcodes = load_barcodes(barcodes_path);
    let barcode_search = BarcodeSearch::new(barcodes, dist_thresh);
    let re_codes: Regex = Regex::new(r".*/1/(.{8})-(.{8})-(.{8})-(.{8})\t").unwrap();
    let lines = open_file(pairs_path).lines();

    // variables shared by threads
    let barcode_search = Arc::new(barcode_search);
    let lines = Arc::new(Mutex::new(lines));
    let re_codes = Arc::new(re_codes);
    let mut counters = vec![];
    for _ in 0..threads { counters.push(Arc::new(Mutex::new(Counter::new()))) };
    let counters = Arc::new(counters);
    let mut handles = vec![];
    let (tx, rx) = mpsc::channel();

    for t_id in 0..threads {
        let lines = Arc::clone(&lines);
        let re_codes = Arc::clone(&re_codes);
        let barcode_search = Arc::clone(&barcode_search);
        let counters = Arc::clone(&counters);
        let tx1 = mpsc::Sender::clone(&tx);
        let handle = thread::spawn(move || {
            let mut counter = Counter::new();
            loop {
                let line = {
                    let mut lines = lines.lock().unwrap();
                    match lines.next() {
                        Some(line) => {
                            line.unwrap()
                        },
                        None => {
                            break
                        }
                    }
                };
                if line.starts_with("#") { continue }
                let rec = PairRec::from_line(&line, &re_codes);
                let b = {
                    let mut counter = counters[t_id as usize].lock().unwrap();
                    locate_barcode(&rec, &barcode_search, max_diff_b1b2, max_diff_r1r2, &mut counter)
                };
                tx1.send((line, b)).unwrap();
            }
        });
        handles.push(handle);
    }

    let mut code_to_file = HashMap::new();

    let open_out_file = |code| {
        let file_name = format!("{}.{}.pairs", output_prefix, code);
        io::BufWriter::new(File::create(file_name).unwrap())
    };

    loop {
        match rx.recv_timeout(Duration::from_millis(wait_t)) {
            Ok((line, b)) => {
                if let Some(code) = b {
                    let f = code_to_file.entry(code.clone()).or_insert(open_out_file(code));
                    write!(*f, "{}\n", line).unwrap();
                }
            },
            _ => {
                info!("End split cell.");
                break;
            }
        }
    }

    for handle in handles {
        handle.join().unwrap();
    }

    let mut counter = counters[0].lock().unwrap().clone();
    for i in 1..counters.len() {
         counter = counter + counters[i].lock().unwrap().clone();
    }
    info!("{}", counter);
    let counter_res_path = format!("{}.count.txt", output_prefix);
    let mut counter_res_file = File::create(counter_res_path).unwrap();
    write!(counter_res_file, "{}", counter).unwrap();

}
