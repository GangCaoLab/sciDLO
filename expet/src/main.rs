use std::fs::File;
use std::io::Read;
use std::thread;
use std::sync::mpsc;
use std::sync::{Mutex, Arc};
use std::time::Duration;
use std::collections::HashMap;
use std::ops;
use std::hash::Hash;
use std::fmt;

extern crate bio;
extern crate clap;
extern crate flate2;
extern crate log;
extern crate simple_logger;

use clap::{Arg, App};
use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment;
use bio::io::{fastq};
use bio::io::fastq::Record;
use log::{info};
use flate2::read::GzDecoder;

fn read_fq(fq_path: String) -> fastq::Reader<Box<dyn Read + Send + Sync>> {
    let fq_file: Box<dyn Read + Send + Sync> = if fq_path.ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(fq_path).unwrap()))
    } else {
        Box::new(File::open(fq_path).unwrap())
    };
    fastq::Reader::new(fq_file)
}

#[derive(Clone)]
struct Counter {
     valid: u64,
     r1_not_match: u64,
     r2_not_match: u64,
     p1_too_short: u64,
     p2_too_short: u64,
     p1_too_long: u64,
     p2_too_long: u64,
     p1_add_base: u64,
     p2_add_base: u64,
     total: u64,
     pet1_len_cnts: HashMap<usize, u64>,
     pet2_len_cnts: HashMap<usize, u64>,
     barcode_cnts: HashMap<String, u64>,
}

impl Counter {
     fn new() -> Self {
          Self {
               valid: 0,
               r1_not_match: 0,
               r2_not_match: 0,
               p1_too_short: 0,
               p2_too_short: 0,
               p1_too_long: 0,
               p2_too_long: 0,
               p1_add_base: 0,
               p2_add_base: 0,
               total: 0,
               pet1_len_cnts: HashMap::new(),
               pet2_len_cnts: HashMap::new(),
               barcode_cnts: HashMap::new(),
          }
     }
}

fn add_hashmap<T: Clone + Eq + Hash>(m1: HashMap<T, u64>, m2: HashMap<T, u64>) -> HashMap<T, u64> {
     let mut m = m1.clone();
     for (k, v) in m2 {
          *m.entry(k).or_insert(v) += v;
     }
     return m
}

impl ops::Add<Counter> for Counter {
     type Output = Counter;

     fn add(self, _rhs: Counter) -> Counter {
          Self {
               valid: self.valid + _rhs.valid,
               r1_not_match: self.r1_not_match + _rhs.r1_not_match,
               r2_not_match: self.r2_not_match + _rhs.r2_not_match,
               p1_too_short: self.p1_too_short + _rhs.p1_too_short,
               p2_too_short: self.p2_too_short + _rhs.p2_too_short,
               p1_too_long: self.p1_too_long + _rhs.p1_too_long,
               p2_too_long: self.p2_too_long + _rhs.p2_too_long,
               p1_add_base: self.p1_add_base + _rhs.p1_add_base,
               p2_add_base: self.p2_add_base + _rhs.p2_add_base,
               total: self.total + _rhs.total,
               pet1_len_cnts: add_hashmap(self.pet1_len_cnts, _rhs.pet1_len_cnts),
               pet2_len_cnts: add_hashmap(self.pet2_len_cnts, _rhs.pet2_len_cnts),
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
               r1_not_match\t{}\t{}\n\
               r2_not_match\t{}\t{}\n\
               p1_too_short\t{}\t{}\n\
               p2_too_short\t{}\t{}\n\
               p1_too_long\t{}\t{}\n\
               p2_too_long\t{}\t{}\n\
               p1_add_base\t{}\t{}\n\
               p2_add_base\t{}\t{}\n\
               total reads\t{}\n\n",
               self.valid, ratio(self.valid),
               self.r1_not_match, ratio(self.r1_not_match),
               self.r2_not_match, ratio(self.r2_not_match),
               self.p1_too_short, ratio(self.p1_too_short),
               self.p2_too_short, ratio(self.p2_too_short),
               self.p1_too_long, ratio(self.p1_too_long),
               self.p2_too_long, ratio(self.p2_too_long),
               self.p1_add_base, ratio(self.p1_add_base),
               self.p2_add_base, ratio(self.p2_add_base),
               self.total,
          );
          msg.push_str("PET1 length distribution:\n");
          let mut keys_pet1: Vec<&usize> = self.pet1_len_cnts.keys().collect();
          keys_pet1.sort();
          for k in keys_pet1 {
               msg.push_str(&format!("{}\t{}\n", k, self.pet1_len_cnts.get(k).unwrap()));
          }
          msg.push_str("\n");
          msg.push_str("PET2 length distribution:\n");
          let mut keys_pet2: Vec<&usize> = self.pet2_len_cnts.keys().collect();
          keys_pet2.sort();
          for k in keys_pet2 {
               msg.push_str(&format!("{}\t{}\n", k, self.pet2_len_cnts.get(k).unwrap()));
          }
          msg.push_str("\n");
          msg.push_str("barcodes counts:\n");
          let mut items_bar_cnts: Vec<(&String, &u64)> = self.barcode_cnts.iter().collect();
          items_bar_cnts.sort_by(|t1, t2| t2.1.cmp(t1.1));
          for (barcode, cnt) in items_bar_cnts {
               msg.push_str(&format!("{}\t{}\n", barcode, cnt));
          }
          write!(f, "{}", msg)
     }
}

struct Extractor {
     linker: Vec<u8>,
     enzyme: Vec<String>,
     _enzyme_half: String,
     score_ratio_thresh: f32,
     n_in_linker: usize,
     min_pet_len: usize,
     max_pet_len: usize,
     pet_cut_len: usize,
     is_extract_barcode: bool,
     _barcode_pos: Vec<(usize, usize)>,
}

fn find_n_blocks(linker: &str) -> Vec<(usize, usize)> {
     let mut start = 0;
     let mut end;
     let mut pos = vec![];
     for i in 1..linker.len() {
          let c = linker.chars().nth(i).unwrap();
          let c_p = linker.chars().nth(i-1).unwrap();
          if (c == 'N') && (c_p != 'N') {
               start = i;
          } else if (c != 'N') && (c_p == 'N') {
               end = i;
               pos.push((start, end));
          }
     }
     pos
}

impl Extractor {

     fn new (linker: &str, enzyme: &str, score_ratio_thresh: f32,
             min_pet_len: usize, max_pet_len: usize, pet_cut_len: usize,
             is_extract_barcode: bool) -> Self {
          let n_in_linker: usize = linker.matches("N").count(); 
          let e_parts: Vec<String> = enzyme.split("^").map(|s| s.to_string()).collect();
          if e_parts.len() != 3 { panic!("Enzyme should contain two cut site.") }
          if min_pet_len >= pet_cut_len || pet_cut_len >= max_pet_len {
               panic!("PET length parameters hould in this relationship: min_pet_len < pet_cut_len < max_pet_len")
          }
          let barcode_pos = find_n_blocks(linker);

          Self {
               linker: linker.as_bytes().to_vec(),
               _enzyme_half: format!("{}{}", e_parts[0], e_parts[1]),
               enzyme: e_parts,
               score_ratio_thresh: score_ratio_thresh,
               n_in_linker: n_in_linker,
               min_pet_len: min_pet_len,
               max_pet_len: max_pet_len,
               pet_cut_len: pet_cut_len,
               is_extract_barcode: is_extract_barcode,
               _barcode_pos: barcode_pos,
          }
     }

     fn extract_pet(&self, rec1: Record, rec2: Record, counter: &mut Counter) -> Result<(Record, Record), ()> {
          let score = |a: u8, b: u8| if (a == b) || (a==78) {1i32} else {-1i32};
          let seq1 = rec1.seq();
          let seq2 = rec2.seq();
          let mut aligner = Aligner::with_capacity(seq1.len(), self.linker.len(), -1, -1, score);
          let aln1 = aligner.semiglobal(&self.linker, seq1);
          let aln2 = aligner.semiglobal(&self.linker, seq2);
          let thresh: f32 = (self.linker.len() - self.n_in_linker) as f32 * self.score_ratio_thresh;
          counter.total += 1;

          if ((aln1.score - self.n_in_linker as i32) as f32) < thresh {
               counter.r1_not_match += 1;
               return Err(())
          }
          let mut pet1 = String::from_utf8(seq1[0..aln1.ystart].to_vec()).unwrap();
          if pet1.ends_with(&self._enzyme_half) {
               counter.p1_add_base += 1;
               pet1.push_str(&self.enzyme[2]);
          }
          if pet1.len() < self.min_pet_len {
               counter.p1_too_short += 1;
               return Err(())
          }
          if pet1.len() > self.max_pet_len {
               counter.p1_too_long += 1;
               pet1 = pet1[(pet1.len() - self.pet_cut_len)..pet1.len()].to_string();
          }

          if ((aln2.score - self.n_in_linker as i32) as f32) < thresh {
               counter.r2_not_match += 1;
               return Err(())
          }
          let mut pet2 = String::from_utf8(seq2[0..aln2.ystart].to_vec()).unwrap();
          if pet2.ends_with(&self._enzyme_half) {
               counter.p2_add_base += 1;
               pet2.push_str(&self.enzyme[2]);
          }
          if pet2.len() < self.min_pet_len {
               counter.p2_too_short += 1;
               return Err(())
          }
          if pet2.len() > self.max_pet_len {
               counter.p2_too_long += 1;
               pet2 = pet2[(pet2.len() - self.pet_cut_len)..pet2.len()].to_string();
          }

          counter.valid += 1;
          *counter.pet1_len_cnts.entry(pet1.len()).or_insert(0) += 1;
          *counter.pet2_len_cnts.entry(pet2.len()).or_insert(0) += 1;

          let mut p1_id = rec1.id().to_string();
          let mut p2_id = rec2.id().to_string();
          if self.is_extract_barcode {
               let barcode = self.extract_barcode(seq1, seq2, &aln1, &aln2);
               p1_id = format!("{}/{}", p1_id, barcode);
               p2_id = format!("{}/{}", p2_id, barcode);
               *counter.barcode_cnts.entry(barcode).or_insert(0) += 1;
          }
          let pet1 = Record::with_attrs(&p1_id, None, pet1.as_bytes(), &rec1.qual()[0..pet1.len()]);
          let pet2 = Record::with_attrs(&p2_id, None, pet2.as_bytes(), &rec2.qual()[0..pet2.len()]);
          return Ok((pet1, pet2))
     }

     fn extract_barcode(&self, seq1: &[u8], seq2: &[u8], aln1: &Alignment, aln2: &Alignment) -> String {
          let mut barcodes = Vec::with_capacity(2*self._barcode_pos.len());
          for (s, e) in &self._barcode_pos {
               let b1 = String::from_utf8(seq1[aln1.ystart+s..aln1.ystart+e].to_vec()).unwrap();
               let b2 = String::from_utf8(seq2[aln2.ystart+s..aln2.ystart+e].to_vec()).unwrap();
               barcodes.push(b1);
               barcodes.push(b2);
          }
          let barcode = barcodes.join("-");
          barcode
     }
}

fn new_writers(prefix: &str, barcode: &Option<String>) -> (fastq::Writer<File>, fastq::Writer<File>) {
     let (pet1_out_path, pet2_out_path) = match barcode {
          Some(b) => (format!("{}_{}.pet1.fq", prefix, b), format!("{}_{}.pet2.fq", prefix, b)),
          None => (format!("{}.pet1.fq", prefix), format!("{}.pet2.fq", prefix))
     };
     let pet1_out_f = File::create(pet1_out_path).unwrap();
     let pet2_out_f = File::create(pet2_out_path).unwrap();
     let writer_pet1 = fastq::Writer::new(pet1_out_f);
     let writer_pet2 = fastq::Writer::new(pet2_out_f);
     (writer_pet1, writer_pet2)
}

fn main() {
     simple_logger::init().unwrap();

     let matches = App::new("Extract and counting seq pairs.")
          .arg(Arg::with_name("fq1")
               .required(true)
               .help("Fastq file of reads 2."))
          .arg(Arg::with_name("fq2")
               .required(true)
               .help("Fastq file of reads 2."))
          .arg(Arg::with_name("linker")
               .short("l")
               .long("linker")
               .required(true)
               .takes_value(true)
               .help("The linker sequence(Not incluede enzyme), \
                     if contain barcode, use 'N' mark the barcode sequence. \
                     like: GTCGGANNNNNNNNGCTAGCNNNNNNNNTCCGAC"))
          .arg(Arg::with_name("split_barcode")
               .short("b")
               .long("split_barcode")
               .takes_value(false)
               .help("Split outputs by barcodes if specify."))
          .arg(Arg::with_name("min_pet_len")
               .short("m")
               .long("min_pet_len")
               .takes_value(true)
               .help("Min length of PET, will be drop if short than this."))
          .arg(Arg::with_name("max_pet_len")
               .short("M")
               .long("max_pet_len")
               .takes_value(true)
               .help("Max length of PET, will be cut if longer than this."))
          .arg(Arg::with_name("pet_cut_len")
               .short("c")
               .long("pet_cut_len")
               .takes_value(true)
               .help("If PET length large than the upper len range, will cut to this length."))
           .arg(Arg::with_name("score_ratio_thresh")
               .short("s")
               .long("score_ratio_thresh")
               .takes_value(true)
               .help("Threshold of (align score / pattern length)"))
          .arg(Arg::with_name("enzyme")
               .short("e")
               .long("enzyme")
               .required(true)
               .takes_value(true)
               .help("Enzyme recognize site, use '^' indicate the cut site, for example T^TA^A"))
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

     let fq1_path = matches.value_of("fq1").unwrap();
     let fq2_path = matches.value_of("fq2").unwrap();
     let linker = matches.value_of("linker").unwrap();
     let split_barcode = matches.is_present("split_barcode");
     let min_pet_len = matches.value_of("min_pet_len").unwrap_or("10");
     let min_pet_len: usize = min_pet_len.parse().unwrap();
     let max_pet_len = matches.value_of("max_pet_len").unwrap_or("22");
     let max_pet_len: usize = max_pet_len.parse().unwrap();
     let pet_cut_len = matches.value_of("max_pet_len").unwrap_or("20");
     let pet_cut_len: usize = pet_cut_len.parse().unwrap();
     let score_ratio_thresh = matches.value_of("score_ratio_thresh").unwrap_or("0.80");
     let score_ratio_thresh: f32 = score_ratio_thresh.parse().unwrap();
     let enzyme = matches.value_of("enzyme").unwrap();
     let output_prefix = matches.value_of("output_prefix").unwrap();
     let threads = matches.value_of("threads").unwrap_or("1");
     let threads: u8 = threads.parse().unwrap();
     let wait_t = matches.value_of("wait_timeout").unwrap_or("500");
     let wait_t: u64 = wait_t.parse().unwrap();

     info!("fastq1: {} fastq2: {}\n linker: {} enzyme: {}\nscore_ratio_thresh: {}, threads: {}",
           fq1_path, fq2_path, linker, enzyme, score_ratio_thresh, threads);

     let recs_1 = read_fq(fq1_path.to_string()).records();
     let recs_2 = read_fq(fq2_path.to_string()).records();

     let extractor = Extractor::new(&linker, &enzyme, score_ratio_thresh,
          min_pet_len, max_pet_len, pet_cut_len, split_barcode);

     let recs_1 = Arc::new(Mutex::new(recs_1));
     let recs_2 = Arc::new(Mutex::new(recs_2));
     let extractor = Arc::new(extractor);
     let mut counters = vec![];
     for _ in 0..threads { counters.push(Arc::new(Mutex::new(Counter::new()))) };
     let counters = Arc::new(counters);
     let mut handles = vec![];
     let (tx, rx) = mpsc::channel();

     for t_id in 0..threads {
          let recs_1 = Arc::clone(&recs_1);
          let recs_2 = Arc::clone(&recs_2);
          let extractor = Arc::clone(&extractor);
          let counters = Arc::clone(&counters);
          let tx1 = mpsc::Sender::clone(&tx);
          let handle = thread::spawn(move || {
               loop {
                    let rec1 = {
                         let mut recs_1 = recs_1.lock().unwrap();
                         match recs_1.next() {
                              Some(r) => match r {
                                  Ok(r_) => r_,
                                  Err(e) => panic!("{:?}", e),
                              },
                              None => break
                         }
                    };
                    let rec2 = {
                         let mut recs_2 = recs_2.lock().unwrap();
                         match recs_2.next() {
                             Some(r) => match r {
                                 Ok(r_) => r_,
                                 Err(e) => panic!("{:?}", e),
                             },
                             None => break
                         }
                    };
                    let res = {
                         let mut counter = counters[t_id as usize].lock().unwrap();
                         extractor.extract_pet(rec1, rec2, &mut counter)
                    };
                    match res {
                         Ok((pet1, pet2)) => {
                              tx1.send((pet1, pet2)).unwrap();
                         },
                         _ => {}
                    }
               }
          });
          handles.push(handle);
     }

     let (mut writer_pet1, mut writer_pet2) = new_writers(output_prefix, &None);
     loop {
          match rx.recv_timeout(Duration::from_millis(wait_t)) {
               Ok((pet1, pet2)) => {
                    writer_pet1.write_record(&pet1).unwrap();
                    writer_pet2.write_record(&pet2).unwrap();
               },
               _ => {
                    info!("End extract PETs.");
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

}
