use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::ops;
use std::sync::mpsc;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

extern crate bio;
extern crate clap;
extern crate log;
extern crate simple_logger;

use bio::alignment::pairwise::Aligner;
use bio::alignment::Alignment;
use bio::alphabets::dna::revcomp;
use bio::io::fastq;
use bio::io::fastq::Record;
use clap::{App, Arg};
use log::info;

use utils::{add_hashmap, open_file};

fn open_fq(fq_path: String) -> fastq::Reader<io::BufReader<Box<dyn Read + Send + Sync>>> {
    let f = open_file(&fq_path);
    fastq::Reader::new(f)
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
    adapter_not_match_rec1: u64,
    adapter_not_match_rec2: u64,
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
            adapter_not_match_rec1: 0,
            adapter_not_match_rec2: 0,
            total: 0,
            pet1_len_cnts: HashMap::new(),
            pet2_len_cnts: HashMap::new(),
            barcode_cnts: HashMap::new(),
        }
    }
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
            adapter_not_match_rec1: self.adapter_not_match_rec1 + _rhs.adapter_not_match_rec1,
            adapter_not_match_rec2: self.adapter_not_match_rec2 + _rhs.adapter_not_match_rec2,
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
            if self.total == 0 {
                return format!("0%");
            }
            format!("{:.2}%", ((c * 100) as f64) / (self.total as f64))
        };
        let mut msg = format!(
            "Count result:\n\
            valid:\t{}\t{}\n\
            r1_not_match\t{}\t{}\n\
            r2_not_match\t{}\t{}\n\
            p1_too_short\t{}\t{}\n\
            p2_too_short\t{}\t{}\n\
            p1_too_long\t{}\t{}\n\
            p2_too_long\t{}\t{}\n\
            p1_add_base\t{}\t{}\n\
            p2_add_base\t{}\t{}\n\
            adapter_not_match_rec1\t{}\t{}\n\
            adapter_not_match_rec2\t{}\t{}\n\
            total reads\t{}\n\n",
            self.valid,
            ratio(self.valid),
            self.r1_not_match,
            ratio(self.r1_not_match),
            self.r2_not_match,
            ratio(self.r2_not_match),
            self.p1_too_short,
            ratio(self.p1_too_short),
            self.p2_too_short,
            ratio(self.p2_too_short),
            self.p1_too_long,
            ratio(self.p1_too_long),
            self.p2_too_long,
            ratio(self.p2_too_long),
            self.p1_add_base,
            ratio(self.p1_add_base),
            self.p2_add_base,
            ratio(self.p2_add_base),
            self.adapter_not_match_rec1,
            ratio(self.adapter_not_match_rec1),
            self.adapter_not_match_rec2,
            ratio(self.adapter_not_match_rec2),
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
    adapter: Option<Vec<u8>>,
    score_ratio_thresh_adapter: f32,
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
        let c_p = linker.chars().nth(i - 1).unwrap();
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
    fn new(
        linker: &str,
        enzyme: &str,
        score_ratio_thresh: f32,
        adapter: Option<&str>,
        score_ratio_thresh_adapter: f32,
        min_pet_len: usize,
        max_pet_len: usize,
        pet_cut_len: usize,
        is_extract_barcode: bool,
    ) -> Self {
        let n_in_linker: usize = linker.matches("N").count();
        let e_parts: Vec<String> = enzyme.split("^").map(|s| s.to_string()).collect();
        if e_parts.len() != 3 {
            panic!("Enzyme should contain two cut site.")
        }
        if min_pet_len >= pet_cut_len || pet_cut_len >= max_pet_len {
            panic!("PET length parameters hould in this relationship: min_pet_len < pet_cut_len < max_pet_len")
        }
        let barcode_pos = find_n_blocks(linker);
        let adapter = match adapter {
            None => None,
            Some(s) => Some(s.as_bytes().to_vec()),
        };

        Self {
            linker: linker.as_bytes().to_vec(),
            _enzyme_half: format!("{}{}", e_parts[0], e_parts[1]),
            enzyme: e_parts,
            score_ratio_thresh: score_ratio_thresh,
            adapter: adapter,
            score_ratio_thresh_adapter: score_ratio_thresh_adapter,
            n_in_linker: n_in_linker,
            min_pet_len: min_pet_len,
            max_pet_len: max_pet_len,
            pet_cut_len: pet_cut_len,
            is_extract_barcode: is_extract_barcode,
            _barcode_pos: barcode_pos,
        }
    }

    fn extract_pet(
        &self,
        rec1: Record,
        rec2: Option<Record>,
        counter: &mut Counter,
    ) -> Result<(Record, Record), ()> {
        counter.total += 1;
        let score = |a: u8, b: u8| if (a == b) || (a == 78) { 1i32 } else { -1i32 };
        let mut seq1 = rec1.seq();

        let mut aligner_ada = None;
        if let Some(adapter) = &self.adapter {
            // Trim adapter in rec1
            let mut aligner_ = Aligner::with_capacity(seq1.len(), adapter.len(), -1, -1, score);
            let aln1_ada = aligner_.semiglobal(&adapter, seq1);
            if aln1_ada.score as f32 >= (self.score_ratio_thresh_adapter * adapter.len() as f32) {
                seq1 = &seq1[0..aln1_ada.ystart];
            } else {
                counter.adapter_not_match_rec1 += 1;
            }
            aligner_ada = Some(aligner_);
        };

        // Align linker to rec1
        let mut aligner = Aligner::with_capacity(seq1.len(), self.linker.len(), -1, -1, score);
        let aln1 = aligner.semiglobal(&self.linker, seq1);
        let thresh: f32 = (self.linker.len() - self.n_in_linker) as f32 * self.score_ratio_thresh;
        if ((aln1.score - self.n_in_linker as i32) as f32) < thresh {
            counter.r1_not_match += 1;
            return Err(());
        }
        // Extract pet1 from rec1's head
        let mut pet1 = String::from_utf8(seq1[0..aln1.ystart].to_vec()).unwrap();
        let mut qual1 = rec1.qual()[0..pet1.len()].to_vec();
        if pet1.ends_with(&self._enzyme_half) {
            // Add adition base to pet1
            counter.p1_add_base += 1;
            pet1.push_str(&self.enzyme[2]);
            qual1.push(70)
        }
        if pet1.len() < self.min_pet_len {
            counter.p1_too_short += 1;
            return Err(());
        }
        if pet1.len() > self.max_pet_len {
            // cut pet1
            counter.p1_too_long += 1;
            pet1 = pet1[(pet1.len() - self.pet_cut_len)..pet1.len()].to_string();
            qual1 = qual1[0..pet1.len()].to_vec()
        }

        let mut pet2;
        let barcode;
        let mut qual2;
        if let Some(rec2) = rec2 {
            // In PE mode, align linker to rec2
            let mut seq2 = rec2.seq();
            if let Some(adapter) = &self.adapter {
                // Trim adapter in rec2
                let mut aligner_ = aligner_ada.unwrap();
                let aln2_ada = aligner_.semiglobal(&adapter, seq2);
                if aln2_ada.score as f32 >= (self.score_ratio_thresh_adapter * adapter.len() as f32)
                {
                    seq2 = &seq2[0..aln2_ada.ystart];
                } else {
                    counter.adapter_not_match_rec2 += 1;
                }
            }
            let aln2 = aligner.semiglobal(&self.linker, seq2);
            if ((aln2.score - self.n_in_linker as i32) as f32) < thresh {
                counter.r2_not_match += 1;
                return Err(());
            }
            // PE mode, extract pet2 from rec2's head
            pet2 = String::from_utf8(seq2[0..aln2.ystart].to_vec()).unwrap();
            qual2 = rec2.qual()[0..pet2.len()].to_vec();
            match self.cut_pet2(&mut pet2, &mut qual2, counter) {
                Ok(()) => {},
                Err(()) => return Err(()),
            };
            barcode = self.extract_barcode_pe(seq1, seq2, &aln1, &aln2);
        } else {
            // SE mode, extract pet2 from rec1's tail
            let p2 = seq1[aln1.yend..seq1.len()].to_vec();
            let p2_rc = revcomp(&p2);
            pet2 = String::from_utf8(p2_rc).unwrap();
            qual2 = rec1.qual()[aln1.yend..seq1.len()].to_vec();
            qual2.reverse();
            match self.cut_pet2(&mut pet2, &mut qual2, counter) {
                Ok(()) => {},
                Err(()) => return Err(()),
            };
            barcode = self.extract_barcode_se(seq1, &aln1);
        }

        counter.valid += 1;
        *counter.pet1_len_cnts.entry(pet1.len()).or_insert(0) += 1;
        *counter.pet2_len_cnts.entry(pet2.len()).or_insert(0) += 1;

        // construct fq records of pet1, pet2
        let mut p_id = rec1.id().to_string();
        if self.is_extract_barcode {
            p_id = format!("{}/{}", p_id, barcode);
            *counter.barcode_cnts.entry(barcode).or_insert(0) += 1;
        }
        let pet1 = Record::with_attrs(&p_id, None, pet1.as_bytes(), &qual1);
        let pet2 = Record::with_attrs(&p_id, None, pet2.as_bytes(), &qual2);
        return Ok((pet1, pet2));
    }

    fn cut_pet2(
        &self,
        pet2: &mut String,
        qual2: &mut Vec<u8>,
        counter: &mut Counter,
    ) -> Result<(), ()> {
        if pet2.ends_with(&self._enzyme_half) {
            // add addition base to pet2
            counter.p2_add_base += 1;
            pet2.push_str(&self.enzyme[2]);
            qual2.push(70);
        }
        if pet2.len() < self.min_pet_len {
            counter.p2_too_short += 1;
            return Err(());
        }
        if pet2.len() > self.max_pet_len {
            // cut pet2
            counter.p2_too_long += 1;
            *pet2 = pet2[(pet2.len() - self.pet_cut_len)..pet2.len()].to_string();
            *qual2 = qual2[0..pet2.len()].to_vec()
        }
        return Ok(());
    }

    fn extract_barcode_pe(
        &self,
        seq1: &[u8],
        seq2: &[u8],
        aln1: &Alignment,
        aln2: &Alignment,
    ) -> String {
        let mut barcodes = Vec::with_capacity(2 * self._barcode_pos.len());
        for (s, e) in &self._barcode_pos {
            let b1 = String::from_utf8(seq1[aln1.ystart + s..aln1.ystart + e].to_vec()).unwrap();
            let b2 = String::from_utf8(seq2[aln2.ystart + s..aln2.ystart + e].to_vec()).unwrap();
            barcodes.push(b1);
            barcodes.push(b2);
        }
        let barcode = barcodes.join("-");
        barcode
    }

    fn extract_barcode_se(&self, seq1: &[u8], aln1: &Alignment) -> String {
        let mut barcodes = Vec::with_capacity(2 * self._barcode_pos.len());
        for (s, e) in &self._barcode_pos {
            let b1 = String::from_utf8(seq1[aln1.ystart + s..aln1.ystart + e].to_vec()).unwrap();
            barcodes.push(b1.clone());
            barcodes.push(b1);
        }
        let barcode = barcodes.join("-");
        barcode
    }
}

fn new_writers(
    prefix: &str,
    barcode: &Option<String>,
) -> (
    fastq::Writer<io::BufWriter<File>>,
    fastq::Writer<io::BufWriter<File>>,
) {
    let (pet1_out_path, pet2_out_path) = match barcode {
        Some(b) => (
            format!("{}_{}.pet1.fq", prefix, b),
            format!("{}_{}.pet2.fq", prefix, b),
        ),
        None => (format!("{}.pet1.fq", prefix), format!("{}.pet2.fq", prefix)),
    };
    let pet1_out_f = io::BufWriter::new(File::create(pet1_out_path).unwrap());
    let pet2_out_f = io::BufWriter::new(File::create(pet2_out_path).unwrap());
    let writer_pet1 = fastq::Writer::new(pet1_out_f);
    let writer_pet2 = fastq::Writer::new(pet2_out_f);
    (writer_pet1, writer_pet2)
}

fn main() {
    simple_logger::SimpleLogger::new().init().unwrap();

    let matches = App::new("Extract PETs.")
        .arg(
            Arg::with_name("fq1")
                .long("fq1")
                .takes_value(true)
                .required(true)
                .help("Fastq file of reads 2."),
        )
        .arg(
            Arg::with_name("fq2")
                .long("fq2")
                .takes_value(true)
                .help("Fastq file of reads 2."),
        )
        .arg(
            Arg::with_name("linker")
                .short("l")
                .long("linker")
                .required(true)
                .takes_value(true)
                .help(
                    "The linker sequence(Not incluede enzyme), \
                     if contain barcode, use 'N' mark the barcode sequence. \
                     like: GTCGGANNNNNNNNGCTAGCNNNNNNNNTCCGAC",
                ),
        )
        .arg(
            Arg::with_name("split_barcode")
                .short("b")
                .long("split_barcode")
                .takes_value(false)
                .help("Split outputs by barcodes if specify."),
        )
        .arg(
            Arg::with_name("min_pet_len")
                .short("m")
                .long("min_pet_len")
                .takes_value(true)
                .help("Min length of PET, will be drop if short than this."),
        )
        .arg(
            Arg::with_name("max_pet_len")
                .short("M")
                .long("max_pet_len")
                .takes_value(true)
                .help("Max length of PET, will be cut if longer than this."),
        )
        .arg(
            Arg::with_name("pet_cut_len")
                .short("c")
                .long("pet_cut_len")
                .takes_value(true)
                .help("If PET length large than the upper len range, will cut to this length."),
        )
        .arg(
            Arg::with_name("score_ratio_thresh")
                .short("s")
                .long("score_ratio_thresh")
                .takes_value(true)
                .help("Threshold of (align score / pattern length)"),
        )
        .arg(
            Arg::with_name("adapter")
                .short("a")
                .long("adapter")
                .takes_value(true)
                .help("The adapter sequence for trim."),
        )
        .arg(
            Arg::with_name("score_ratio_thresh_adapter")
                .long("score_ratio_thresh_adapter")
                .takes_value(true)
                .help("Threshold of adapter's (align score / pattern length)"),
        )
        .arg(
            Arg::with_name("enzyme")
                .short("e")
                .long("enzyme")
                .required(true)
                .takes_value(true)
                .help("Enzyme recognize site, use '^' indicate the cut site, for example T^TA^A"),
        )
        .arg(
            Arg::with_name("output_prefix")
                .short("o")
                .long("output_prefix")
                .required(true)
                .takes_value(true)
                .help("Prefix of output files."),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .help("Number of threads used for processing reads."),
        )
        .arg(
            Arg::with_name("wait_timeout")
                .long("wait_timeout")
                .takes_value(true)
                .help("Wait time for end channel timeout."),
        )
        .get_matches();

    let fq1_path = matches.value_of("fq1").unwrap();
    let pe_mode = matches.is_present("fq2");
    let fq2_path = if pe_mode {
        Some(matches.value_of("fq2").unwrap())
    } else {
        None
    };
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
    let adapter = if matches.is_present("adapter") {
        Some(matches.value_of("adapter").unwrap())
    } else {
        None
    };
    let sr_th_adapter = matches
        .value_of("score_ratio_thresh_adapter")
        .unwrap_or("0.80");
    let sr_th_adapter: f32 = sr_th_adapter.parse().unwrap();
    let enzyme = matches.value_of("enzyme").unwrap();
    let output_prefix = matches.value_of("output_prefix").unwrap();
    let threads = matches.value_of("threads").unwrap_or("1");
    let threads: u8 = threads.parse().unwrap();
    let wait_t = matches.value_of("wait_timeout").unwrap_or("2000");
    let wait_t: u64 = wait_t.parse().unwrap();

    info!(
        "fastq1: {} fastq2: {:?} pe_mode: {}\n\
          linker: {} enzyme: {} score_ratio_thresh: {}\n\
          adapter: {:?} score_ratio_thresh_adapter: {}\n\
          threads: {}",
        fq1_path,
        fq2_path,
        pe_mode,
        linker,
        enzyme,
        score_ratio_thresh,
        adapter,
        sr_th_adapter,
        threads
    );

    let recs_1 = open_fq(fq1_path.to_string()).records();
    let recs_2 = match fq2_path {
        None => None,
        Some(path) => Some(open_fq(path.to_string()).records()),
    };

    let extractor = Extractor::new(
        &linker,
        &enzyme,
        score_ratio_thresh,
        adapter,
        sr_th_adapter,
        min_pet_len,
        max_pet_len,
        pet_cut_len,
        split_barcode,
    );

    let recs = Arc::new(Mutex::new((recs_1, recs_2)));
    let extractor = Arc::new(extractor);
    let mut counters = vec![];
    for _ in 0..threads {
        counters.push(Arc::new(Mutex::new(Counter::new())))
    }
    let counters = Arc::new(counters);
    let mut handles = vec![];
    let (tx, rx) = mpsc::channel();

    for t_id in 0..threads {
        let recs = Arc::clone(&recs);
        let extractor = Arc::clone(&extractor);
        let counters = Arc::clone(&counters);
        let tx1 = mpsc::Sender::clone(&tx);
        let handle = thread::spawn(move || loop {
            let (rec1, rec2) = {
                let mut recs = recs.lock().unwrap();
                let rec1 = match recs.0.next() {
                    Some(r) => match r {
                        Ok(r_) => r_,
                        Err(e) => panic!("{:?}", e),
                    },
                    None => break,
                };
                let rec2 = match &mut recs.1 {
                    Some(rec) => match rec.next() {
                        Some(r) => match r {
                            Ok(r_) => Some(r_),
                            Err(e) => panic!("{:?}", e),
                        },
                        None => break,
                    },
                    None => None,
                };
                (rec1, rec2)
            };
            let res = {
                let mut counter = counters[t_id as usize].lock().unwrap();
                extractor.extract_pet(rec1, rec2, &mut counter)
            };
            match res {
                Ok((pet1, pet2)) => {
                    tx1.send((pet1, pet2)).unwrap();
                }
                _ => {}
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
            }
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

    let counter_res_path = format!("{}.count.txt", output_prefix);
    let mut counter_res_file = File::create(counter_res_path).unwrap();
    write!(counter_res_file, "{}", counter).unwrap();
}
