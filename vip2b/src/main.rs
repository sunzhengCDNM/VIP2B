/// vip2b — unified 2bRAD virome profiler
///
/// Rust reimplementation of VIP2B.py (digest + abundance in one command).
/// Parameter names and defaults match the original Python orchestrator exactly.
///
/// Sample list (-i) format — same as original VIP2B.py:
///   Single-end (2 cols):  sample_id <TAB> reads.fastq.gz
///   Paired-end (3 cols):  sample_id <TAB> reads_R1.fastq.gz <TAB> reads_R2.fastq.gz
///   Lines beginning with # are ignored.
///
/// Species identification threshold (-t):
///   G{n}   — G-score > n: implemented in Rust (e.g. G2, G5)
///   M{p}   — ML probability: MAP2B_ML.py called automatically as Phase 3 subprocess.
///
/// The marisa database (.marisa) is read directly by Python — no conversion needed.
/// -c (cutoff) and --intersection apply to mkdb (CreatDB4AllLevel.py) only.

use std::{
    env,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    sync::{
        atomic::{AtomicU64, Ordering::Relaxed},
        Mutex,
    },
};

use clap::Parser;
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use rayon::prelude::*;

// ═══════════════════════════════════════════════════════════════════════════════
// §1  Nucleotide tables
// ═══════════════════════════════════════════════════════════════════════════════

const A_BIT: u8 = 0b0001;
const C_BIT: u8 = 0b0010;
const G_BIT: u8 = 0b0100;
const T_BIT: u8 = 0b1000;
const ANY: u8   = 0b1111;

fn build_nuc_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = A_BIT; t[b'a' as usize] = A_BIT;
    t[b'C' as usize] = C_BIT; t[b'c' as usize] = C_BIT;
    t[b'G' as usize] = G_BIT; t[b'g' as usize] = G_BIT;
    t[b'T' as usize] = T_BIT; t[b't' as usize] = T_BIT;
    t
}

fn build_comp_table() -> [u8; 256] {
    let mut t = [b'N'; 256];
    for (f, to) in [
        (b'A', b'T'), (b'a', b'T'), (b'T', b'A'), (b't', b'A'),
        (b'G', b'C'), (b'g', b'C'), (b'C', b'G'), (b'c', b'G'),
        (b'N', b'N'), (b'n', b'N'),
    ] { t[f as usize] = to; }
    t
}

#[inline]
fn rev_comp(seq: &[u8], comp: &[u8; 256]) -> Vec<u8> {
    seq.iter().rev().map(|&c| comp[c as usize]).collect()
}

#[inline]
fn to_upper(c: u8) -> u8 {
    if c >= b'a' && c <= b'z' { c - 32 } else { c }
}

// ═══════════════════════════════════════════════════════════════════════════════
// §2  Enzyme pattern matching
// ═══════════════════════════════════════════════════════════════════════════════

struct Enzyme {
    name:   &'static str,
    mask:   Vec<u8>,
    len:    usize,
    anchor: Option<(usize, u8)>,
}

impl Enzyme {
    fn new(name: &'static str, regex: &'static str) -> Self {
        let mask = parse_enzyme_regex(regex);
        let len  = mask.len();
        let anchor = mask.iter().enumerate()
            .filter(|(_, &m)| m != ANY)
            .min_by_key(|(_, &m)| m.count_ones())
            .map(|(i, &m)| (i, m));
        Enzyme { name, mask, len, anchor }
    }

    #[inline]
    fn matches(&self, seq: &[u8], nuc: &[u8; 256]) -> bool {
        if let Some((apos, amask)) = self.anchor {
            if (nuc[seq[apos] as usize] & amask) == 0 { return false; }
        }
        self.mask.iter().zip(&seq[..self.len])
            .all(|(&m, &c)| (nuc[c as usize] & m) != 0)
    }
}

fn parse_enzyme_regex(regex: &'static str) -> Vec<u8> {
    let inner = if regex.starts_with("(?=(") && regex.ends_with("))") {
        &regex[4..regex.len() - 2]
    } else { regex };

    let b = inner.as_bytes();
    let mut mask = Vec::with_capacity(40);
    let mut i = 0;
    while i < b.len() {
        let char_mask: u8;
        if b[i] == b'[' {
            i += 1;
            let mut m = 0u8;
            while i < b.len() && b[i] != b']' {
                m |= match b[i] { b'A' => A_BIT, b'C' => C_BIT, b'G' => G_BIT, b'T' => T_BIT, _ => 0 };
                i += 1;
            }
            i += 1;
            char_mask = m;
        } else if b[i].is_ascii_alphabetic() {
            char_mask = match b[i] { b'A' => A_BIT, b'C' => C_BIT, b'G' => G_BIT, b'T' => T_BIT, _ => 0 };
            i += 1;
        } else { i += 1; continue; }

        let repeat: usize = if i < b.len() && b[i] == b'{' {
            i += 1;
            let start = i;
            while i < b.len() && b[i] != b'}' { i += 1; }
            let n: usize = std::str::from_utf8(&b[start..i]).unwrap().parse().unwrap();
            i += 1;
            n
        } else { 1 };
        for _ in 0..repeat { mask.push(char_mask); }
    }
    mask
}

const ALL_ENZYME_DEFS: &[(&str, &str)] = &[
    ("AlfI",   r"(?=([AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}))"),
    ("AloI",   r"(?=([AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}))"),
    ("BaeI",   r"(?=([AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}))"),
    ("BcgI",   r"(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))"),
    ("BplI",   r"(?=([AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}))"),
    ("BsaXI",  r"(?=([AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}))"),
    ("BslFI",  r"(?=([AGCT]{6}GGGAC[AGCT]{14}))"),
    ("Bsp24I", r"(?=([AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}))"),
    ("CjeI",   r"(?=([AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}))"),
    ("CjePI",  r"(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))"),
    ("CspCI",  r"(?=([AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}))"),
    ("FalI",   r"(?=([AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}))"),
    ("HaeIV",  r"(?=([AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}))"),
    ("Hin4I",  r"(?=([AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}))"),
    ("PpiI",   r"(?=([AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}))"),
    ("PsrI",   r"(?=([AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}))"),
];

/// The standard 8-enzyme set that matches the default database.
/// Used when -e all is specified (same as the default).
const DEFAULT_8_ENZYMES: &[&str] = &[
    "AlfI", "BcgI", "BslFI", "CjeI", "CjePI", "FalI", "HaeIV", "Hin4I",
];

fn select_enzymes(enzyme_arg: &str) -> Vec<Enzyme> {
    let wanted: Vec<&str> = if enzyme_arg.trim() == "all" {
        DEFAULT_8_ENZYMES.to_vec()
    } else {
        enzyme_arg.split(',').map(str::trim).collect()
    };
    ALL_ENZYME_DEFS.iter()
        .filter(|(name, _)| wanted.contains(name))
        .map(|(name, regex)| Enzyme::new(name, regex))
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════════
// §3  Streaming FASTQ / FASTA reader
// ══════════════════════════════════════════════════════════════════════════════

fn open_reader(path: &str) -> Box<dyn BufRead + Send> {
    let f = File::open(path).unwrap_or_else(|e| panic!("Cannot open {path}: {e}"));
    if path.ends_with(".gz") {
        Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, f))
    }
}

struct RecordReader {
    reader:     Box<dyn BufRead + Send>,
    is_fastq:   bool,
    line:       String,
    pending_id: Option<String>,
    eof:        bool,
}

impl RecordReader {
    fn new(mut reader: Box<dyn BufRead + Send>) -> Self {
        let mut line = String::with_capacity(256);
        let eof = reader.read_line(&mut line).unwrap_or(0) == 0;
        let is_fastq = line.starts_with('@');
        let pending_id = if eof { None } else { Some(header_id(&line)) };
        RecordReader { reader, is_fastq, line, pending_id, eof }
    }

    fn next_record(&mut self) -> Option<(String, Vec<u8>)> {
        if self.eof { return None; }
        if self.is_fastq {
            let id = self.pending_id.take()?;
            self.line.clear();
            if self.reader.read_line(&mut self.line).unwrap_or(0) == 0 { self.eof = true; return None; }
            let seq: Vec<u8> = self.line.trim_end().bytes().map(to_upper).collect();
            self.line.clear(); self.reader.read_line(&mut self.line).unwrap_or(0);
            self.line.clear(); self.reader.read_line(&mut self.line).unwrap_or(0);
            self.line.clear();
            if self.reader.read_line(&mut self.line).unwrap_or(0) == 0 { self.eof = true; }
            else { self.pending_id = Some(header_id(&self.line)); }
            Some((id, seq))
        } else {
            let id = self.pending_id.take()?;
            let mut seq = Vec::with_capacity(200);
            loop {
                self.line.clear();
                if self.reader.read_line(&mut self.line).unwrap_or(0) == 0 { self.eof = true; break; }
                if self.line.starts_with('>') { self.pending_id = Some(header_id(&self.line)); break; }
                seq.extend(self.line.trim_end().bytes().map(to_upper));
            }
            Some((id, seq))
        }
    }

    fn next_batch(&mut self, n: usize) -> Vec<(String, Vec<u8>)> {
        let mut batch = Vec::with_capacity(n);
        while batch.len() < n {
            match self.next_record() { Some(r) => batch.push(r), None => break }
        }
        batch
    }
}

fn header_id(line: &str) -> String {
    line.trim_start_matches(|c| c == '@' || c == '>')
        .split_whitespace().next().unwrap_or("").to_string()
}

// ═══════════════════════════════════════════════════════════════════════════════
// §4  Digest logic
// ══════════════════════════════════════════════════════════════════════════════

struct TagRecord { id: String, seq: Vec<u8> }

fn format_tag_id(read_id: &str, strand: char, pos: usize, enzyme_name: &str, value_len: usize) -> String {
    let full = format!("{read_id}_{strand}_{pos}_{enzyme_name}");
    let justified = if full.len() >= value_len {
        full[full.len() - value_len..].to_string()
    } else {
        format!("{}{full}", ".".repeat(value_len - full.len()))
    };
    justified.trim_start_matches('.').to_string()
}

fn process_read(
    read_id: &str, seq: &[u8],
    enzymes: &[Enzyme], nuc: &[u8; 256], comp: &[u8; 256], value_len: usize,
) -> (Vec<TagRecord>, u64, u64) {
    let rc = rev_comp(seq, comp);
    let mut records = Vec::new();
    let mut ori_count = 0u64;
    let mut dig_count = 0u64;
    for enzyme in enzymes {
        ori_count += 1;
        let plen = enzyme.len;
        if seq.len() >= plen {
            for i in 0..=(seq.len() - plen) {
                if enzyme.matches(&seq[i..], nuc) {
                    dig_count += 1;
                    records.push(TagRecord {
                        seq: seq[i..i+plen].to_vec(),
                        id:  format_tag_id(read_id, '+', i+1, enzyme.name, value_len),
                    });
                }
            }
        }
        if rc.len() >= plen {
            for i in 0..=(rc.len() - plen) {
                if enzyme.matches(&rc[i..], nuc) {
                    dig_count += 1;
                    records.push(TagRecord {
                        seq: rc[i..i+plen].to_vec(),
                        id:  format_tag_id(read_id, '-', i+1, enzyme.name, value_len),
                    });
                }
            }
        }
    }
    (records, ori_count, dig_count)
}

fn digest_file(
    path: &str, enzymes: &[Enzyme], nuc: &[u8; 256], comp: &[u8; 256],
    value_len: usize, batch_size: usize,
    writer: &Mutex<GzEncoder<BufWriter<File>>>,
    total_ori: &AtomicU64, total_dig: &AtomicU64,
) {
    let mut rr = RecordReader::new(open_reader(path));
    loop {
        let batch = rr.next_batch(batch_size);
        if batch.is_empty() { break; }
        let results: Vec<_> = batch.par_iter()
            .map(|(id, seq)| process_read(id, seq, enzymes, nuc, comp, value_len))
            .collect();
        let mut w = writer.lock().unwrap();
        for (records, ori, dig) in results {
            total_ori.fetch_add(ori, Relaxed);
            total_dig.fetch_add(dig, Relaxed);
            for rec in records {
                w.write_all(b">").unwrap();
                w.write_all(rec.id.as_bytes()).unwrap();
                w.write_all(b"\n").unwrap();
                w.write_all(&rec.seq).unwrap();
                w.write_all(b"\n").unwrap();
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// §5  Species identification threshold (-t)
// ═══════════════════════════════════════════════════════════════════════════════

#[derive(Clone)]
enum FilterMode {
    /// G-score threshold applied inside Rust; pred.result written by gscore_filter.py
    GScore(f64),
    /// ML probability: .xls written in full, then MAP2B_ML.py called as subprocess
    ML(f64),
}

fn parse_threshold(s: &str) -> FilterMode {
    if let Some(v) = s.strip_prefix('G').or_else(|| s.strip_prefix('g')) {
        FilterMode::GScore(v.parse().unwrap_or_else(|_| {
            eprintln!("ERROR: -t {s}: expected number after 'G' (e.g. G5, G2)");
            std::process::exit(1);
        }))
    } else if let Some(v) = s.strip_prefix('M').or_else(|| s.strip_prefix('m')) {
        FilterMode::ML(v.parse().unwrap_or_else(|_| {
            eprintln!("ERROR: -t {s}: expected number after 'M' (e.g. M0.5, M0.1)");
            std::process::exit(1);
        }))
    } else {
        eprintln!("ERROR: -t {s}: must start with G or M (e.g. G5 or M0.5)");
        std::process::exit(1);
    }
}

/// Locate a companion directory (e.g. "scripts", "config") by searching
/// several candidate paths relative to the binary and the working directory.
///
/// Search order:
///   1. $VIP2B_HOME/{rel}                  (explicit override)
///   2. {binary}/../../{rel}               (cargo dev build: target/release/vip2b)
///   3. {binary}/../{rel}                  (installed: bin/vip2b)
///   4. ./{rel}                            (cwd fallback)
///
/// Panics with a clear message if nothing is found.
fn find_sibling(rel: &str) -> std::path::PathBuf {
    // 1. Explicit override via environment variable
    if let Ok(home) = env::var("VIP2B_HOME") {
        let p = std::path::PathBuf::from(&home).join(rel);
        if p.exists() { return p; }
    }

    let bin = env::current_exe().unwrap_or_default();
    let bin_dir = bin.parent().unwrap_or(std::path::Path::new("."));

    // 2. Two levels up (target/release/vip2b → project root)
    let p = bin_dir.join("../..").join(rel);
    if p.exists() { return p; }

    // 3. One level up (installed as bin/vip2b)
    let p = bin_dir.join("..").join(rel);
    if p.exists() { return p; }

    // 4. Current working directory
    let p = std::path::PathBuf::from(rel);
    if p.exists() { return p; }

    eprintln!("ERROR: Cannot find '{rel}' next to the binary or in the current directory.");
    eprintln!("       Set VIP2B_HOME to the vip2b-rs project root, e.g.:");
    eprintln!("         export VIP2B_HOME=/path/to/vip2b-rs");
    std::process::exit(1);
}

/// After abundance phase: call MAP2B_ML.py as a subprocess (for -t M).
fn run_map2b_ml(sample_id: &str, xls_path: &str, outdir: &str, threshold: f64) {
    let script = find_sibling("scripts/MAP2B_ML.py");
    let model  = find_sibling("config/XGB_none_0238.pkl");
    let pred   = format!("{outdir}/pred.result");

    eprintln!("INFO [{sample_id}]: MAP2B_ML.py  threshold={threshold}");
    let status = std::process::Command::new("python3")
        .args([
            script.to_str().unwrap_or("scripts/MAP2B_ML.py"),
            "-i", xls_path,
            "-m", model.to_str().unwrap_or("config/XGB_none_0238.pkl"),
            "-n", "none",
            "-T", &threshold.to_string(),
            "-o", &pred,
        ])
        .status();
    match status {
        Ok(s) if s.success() => eprintln!("INFO [{sample_id}]: pred.result → {pred}"),
        Ok(s) => eprintln!("WARN [{sample_id}]: MAP2B_ML.py exited {s}"),
        Err(e) => eprintln!("WARN [{sample_id}]: Cannot run MAP2B_ML.py: {e}\n            HINT: ensure python3, xgboost and scikit-learn are in PATH"),
    }
}

/// After abundance phase: call gscore_filter.py as a subprocess (for -t G).
fn run_gscore_filter(sample_id: &str, xls_path: &str, outdir: &str, threshold: f64) {
    let script = find_sibling("scripts/gscore_filter.py");
    let pred   = format!("{outdir}/pred.result");

    eprintln!("INFO [{sample_id}]: gscore_filter.py  threshold={threshold}");
    let status = std::process::Command::new("python3")
        .args([
            script.to_str().unwrap_or("scripts/gscore_filter.py"),
            "-i", xls_path,
            "-o", &pred,
            "-g", &(threshold as i64).to_string(),
        ])
        .status();
    match status {
        Ok(s) if s.success() => eprintln!("INFO [{sample_id}]: pred.result → {pred}"),
        Ok(s) => eprintln!("WARN [{sample_id}]: gscore_filter.py exited {s}"),
        Err(e) => eprintln!("WARN [{sample_id}]: Cannot run gscore_filter.py: {e}"),
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// §6  Per-sample abundance via Python subprocess
// ══════════════════════════════════════════════════════════════════════════════

fn run_abundance_script(
    parent_dir:  &str,
    sample_id:   &str,
    fa_path:     &str,
    database:    &str,
    classify:    &str,
    level:       &str,
    ct:          f64,
    scripts_dir: &std::path::Path,
) {
    fs::create_dir_all(&format!("{parent_dir}/{sample_id}"))
        .unwrap_or_else(|e| panic!("Cannot create {parent_dir}/{sample_id}: {e}"));

    let list_path = format!("{parent_dir}/{sample_id}/{sample_id}.ct{ct}.list");
    {
        let mut f = BufWriter::new(
            File::create(&list_path).expect("Cannot create sample list"),
        );
        writeln!(f, "{sample_id}\t{fa_path}").unwrap();
    }

    let script = scripts_dir.join("CalculateRelativeAbundance_Single2bEnzyme.py");
    let status = std::process::Command::new("python3")
        .args([
            script.to_str().unwrap_or("scripts/CalculateRelativeAbundance_Single2bEnzyme.py"),
            "-d", database,
            "-c", classify,
            "-l", &list_path,
            "-t", level,
            "-o", parent_dir,
            "-p", "1",
            "-ct", &ct.to_string(),
        ])
        .status();

    match status {
        Ok(s) if s.success() => {},
        Ok(s) => eprintln!("WARN [{sample_id}]: abundance script (ct={ct}) exited {s}"),
        Err(e) => eprintln!("WARN [{sample_id}]: cannot run abundance script: {e}"),
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// §7  Sample list parsing  (matches VIP2B.py -i format)
// ══════════════════════════════════════════════════════════════════════════════

struct SampleEntry { id: String, r1: String, r2: Option<String> }

fn parse_sample_list(path: &str) -> Vec<SampleEntry> {
    let mut samples = Vec::new();
    for (ln, line) in open_reader(path).lines().enumerate() {
        let line = line.expect("I/O error reading sample list");
        let line = line.trim().to_string();
        if line.is_empty() || line.starts_with('#') { continue; }
        let cols: Vec<&str> = line.splitn(3, '\t').collect();
        match cols.len() {
            2 => samples.push(SampleEntry {
                id: cols[0].to_string(), r1: cols[1].to_string(), r2: None }),
            3 => samples.push(SampleEntry {
                id: cols[0].to_string(), r1: cols[1].to_string(),
                r2: Some(cols[2].to_string()) }),
            _ => eprintln!("WARN: line {}: expected 2 or 3 columns — skipped", ln+1),
        }
    }
    samples
}

// ══════════════════════════════════════════════════════════════════════════════
// §8  CLI  — exact parameter parity with VIP2B.py
//
// NOTE: -o and -d use plain relative placeholder strings as default_value so
// that `--help` never exposes the user's absolute directory paths.
// The actual paths are resolved at runtime in main() via resolve_output() and
// resolve_database() below.
// ═══════════════════════════════════════════════════════════════════════════════

/// Sentinel value used in Args for -o when no path is provided.
const DEFAULT_OUTPUT_SENTINEL: &str = "VIP2B_result";

/// Sentinel value used in Args for -d when no path is provided.
const DEFAULT_DATABASE_SENTINEL: &str = "database/8Enzyme";

/// Resolve the output directory: if the user left it at the sentinel, make it
/// absolute relative to the current working directory.
fn resolve_output(raw: &str) -> String {
    if raw == DEFAULT_OUTPUT_SENTINEL {
        format!("{}/VIP2B_result", env::current_dir().unwrap().display())
    } else {
        raw.to_string()
    }
}

/// Resolve the database prefix: if the user left it at the sentinel, locate it
/// relative to the binary (matching the original default_database() behaviour).
fn resolve_database(raw: &str) -> String {
    if raw == DEFAULT_DATABASE_SENTINEL {
        let bin = env::current_exe().unwrap_or_default();
        let bin_dir = bin.parent().unwrap_or(std::path::Path::new("."));
        bin_dir.join("../database/8Enzyme")
            .to_string_lossy().into_owned()
    } else {
        raw.to_string()
    }
}

#[derive(Parser)]
#[command(
    name = "vip2b",
    about = "Unified 2bRAD virome profiler (digest + abundance)",
    long_about = "\
Runs the complete VIP2B pipeline in two internal phases.\n\
Parameters match the original VIP2B.py exactly.\n\
\n\
Sample list (-i) — lines beginning with # are ignored:\n\
  Single-end (2 cols):  sample_id\\treads.fastq.gz\n\
  Paired-end (3 cols):  sample_id\\treads_R1.fastq.gz\\treads_R2.fastq.gz\n\
\n\
Species identification threshold (-t):\n\
  G{n}  G-score > n — applied in Rust             (common: G2, G5)\n\
  M{p}  ML probability — MAP2B_ML.py called automatically (common: M0.1, M0.5)\n\
\n\
Flags -c and --intersection apply to the mkdb step (CreatDB4AllLevel.py)\n\
and are stored here for CLI compatibility only.",
    version = "1.0.0"
)]
struct Args {
    /// Sample list TSV (2-col SE or 3-col PE; # lines ignored)
    /// Format: sample_id <TAB> reads.fastq.gz [<TAB> reads_R2.fastq.gz]
    #[arg(short = 'i', value_name = "FILE")]
    input: String,

    /// Output directory [default: VIP2B_result in current directory]
    #[arg(short = 'o', default_value = DEFAULT_OUTPUT_SENTINEL, value_name = "DIR")]
    output: String,

    /// Taxonomy level [Class|Order|Family|Genus|Species]
    #[arg(short = 'l',
          value_parser = ["Class", "Order", "Family", "Genus", "Species"],
          default_value = "Species",
          value_name = "LEVEL")]
    level: String,

    /// Enzyme(s): comma-separated names or "all".
    /// Available: AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I
    /// The default 8-enzyme set matches the standard database.
    #[arg(short = 'e',
          default_value = "AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I",
          value_name = "ENZYME[,...]")]
    enzyme: String,

    /// Database prefix (without extension, e.g. /path/to/database/8Enzyme)
    /// [default: database/8Enzyme relative to the vip2b binary]
    #[arg(short = 'd', default_value = DEFAULT_DATABASE_SENTINEL, value_name = "PREFIX")]
    database: String,

    /// Number of parallel processes/threads (more may require more memory)
    #[arg(short = 'p', default_value_t = 1, value_name = "N")]
    processes: usize,

    /// Threshold for species identification.
    /// G{n}: G-score filter applied in Rust (e.g. G2, G5).
    /// M{p}: MAP2B_ML.py is called automatically after the abundance step.
    #[arg(short = 't', default_value = "M0.5", value_name = "G5|M0.5")]
    threshold: String,

    /// Cut-off for database building (mkdb step only; not applied in Rust)
    #[arg(short = 'c', default_value_t = 30000, value_name = "N")]
    cutoff: u32,

    /// Threshold for coverage filtering: sequenced_tags / theoretical_tags
    #[arg(short = 'f', default_value_t = 0.6, value_name = "FLOAT")]
    cov_thresh: f64,

    /// Use intersection of tags between genomes (default: union).
    /// Applies to mkdb step only; stored here for pipeline compatibility.
    #[arg(long)]
    intersection: bool,

    // ── Rust-specific performance option (no Python equivalent) ───────────────
    /// Reads per parallel batch in digest phase
    #[arg(long, default_value_t = 100_000, value_name = "N")]
    batch_size: usize,
}

// ══════════════════════════════════════════════════════════════════════════════
// §9  Entry point
// ═══════════════════════════════════════════════════════════════════════════════

fn main() {
    let args = Args::parse();

    // ── Resolve -o and -d to real absolute paths ───────────────────────────────
    // (args.output / args.database may hold sentinel placeholders when the user
    //  did not supply those flags; resolve them here before any path is used.)
    let output   = resolve_output(&args.output).trim_end_matches('/').to_string();
    let database = resolve_database(&args.database);

    // ── Parse -t threshold ────────────────────────────────────────────────────
    let filter_mode = parse_threshold(&args.threshold);
    match &filter_mode {
        FilterMode::GScore(n) =>
            eprintln!("INFO: Species filter: G-score > {n} → gscore_filter.py (Phase 3)"),
        FilterMode::ML(p) =>
            eprintln!("INFO: Species filter: ML probability > {p} → MAP2B_ML.py (Phase 3)"),
    }
    if args.intersection {
        eprintln!("INFO: --intersection noted (mkdb step only; not applied in Rust digest/abundance)");
    }

    // ── Thread pool ────────────────────────────────────────────────────────────
    if args.processes > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.processes)
            .build_global()
            .expect("Failed to build Rayon thread pool");
    }

    // ── Sample list ────────────────────────────────────────────────────────────
    let samples = parse_sample_list(&args.input);
    if samples.is_empty() {
        eprintln!("ERROR: No samples found in {}", args.input); std::process::exit(1);
    }
    let pe_count = samples.iter().filter(|s| s.r2.is_some()).count();
    let se_count = samples.len() - pe_count;
    eprintln!("INFO: {} sample(s) ({se_count} single-end, {pe_count} paired-end)", samples.len());

    // ── Enzyme patterns ────────────────────────────────────────────────────────
    let enzymes = select_enzymes(&args.enzyme);
    if enzymes.is_empty() {
        eprintln!("ERROR: No valid enzymes matched '{}'", args.enzyme); std::process::exit(1);
    }
    eprintln!("INFO: {} enzyme(s): {}",
        enzymes.len(), enzymes.iter().map(|e| e.name).collect::<Vec<_>>().join(", "));

    let nuc_table  = build_nuc_table();
    let comp_table = build_comp_table();

    let digest_dir    = format!("{output}/0.dige");
    let abundance_dir = format!("{output}/1.qual");
    fs::create_dir_all(&digest_dir)   .expect("Cannot create digest dir");
    fs::create_dir_all(&abundance_dir).expect("Cannot create abundance dir");

    // ══════════════════════════════════════════════════════════════════════════
    // PHASE 1 — Tag extraction (digest)
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 1: Tag extraction ══════════════════════════════════════");

    let digest_outputs: Vec<(String, String)> = samples.par_iter().map(|entry| {
        let out_fa    = format!("{digest_dir}/{}.fa.gz",              entry.id);
        let smp_dir   = format!("{digest_dir}/{}", entry.id);
        let out_stat  = format!("{smp_dir}/{}.dige.stat.xls",          entry.id);
        fs::create_dir_all(&smp_dir).expect("Cannot create per-sample digest dir");

        let fa_file = File::create(&out_fa)
            .unwrap_or_else(|e| panic!("Cannot create {out_fa}: {e}"));
        let gz = GzEncoder::new(BufWriter::with_capacity(1 << 20, fa_file), Compression::default());
        let writer    = Mutex::new(gz);
        let total_ori = AtomicU64::new(0);
        let total_dig = AtomicU64::new(0);

        eprintln!("INFO [{smp}]: R1={r1}{pe}", smp=entry.id, r1=entry.r1,
            pe=entry.r2.as_deref().map(|r| format!("  R2={r}")).unwrap_or_default());

        digest_file(&entry.r1, &enzymes, &nuc_table, &comp_table,
                    50, args.batch_size, &writer, &total_ori, &total_dig);
        if let Some(r2) = &entry.r2 {
            digest_file(r2, &enzymes, &nuc_table, &comp_table,
                        50, args.batch_size, &writer, &total_ori, &total_dig);
        }
        writer.into_inner().unwrap().finish().expect("Failed to finalise gzip");

        let ori = total_ori.load(Relaxed);
        let dig = total_dig.load(Relaxed);
        let pct = if ori > 0 { dig as f64 / ori as f64 * 100.0 } else { 0.0 };
        eprintln!("INFO [{smp}]: input={ori}  extracted={dig}  rate={pct:.2}%", smp=entry.id);

        let mut stat = BufWriter::new(File::create(&out_stat).expect("Cannot create stat file"));
        writeln!(stat, "sample\tenzyme\tinput_sequence_num\tenzyme_reads_num\tpercent").unwrap();
        writeln!(stat, "{}\t{}\t{ori}\t{dig}\t{pct:.2}%", entry.id, args.enzyme).unwrap();

        (entry.id.clone(), out_fa)
    }).collect();

    // ═════════════════════════════════════════════════════════════════════════
    // PHASE 2 — Qualitative abundance
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 2: Qualitative abundance ════════════════════════════════");

    let db_parent = std::path::Path::new(&database)
        .parent()
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|| ".".to_string());
    let classify_path = format!("{db_parent}/abfh_classify_with_speciename.txt.gz");
    let scripts_dir   = find_sibling("scripts");

    let qual_db = format!("{}.{}.uniq", database, args.level);
    let qual_dir = format!("{output}/1.qual");
    fs::create_dir_all(&qual_dir).expect("Cannot create qual dir");

    digest_outputs.par_iter().for_each(|(sample_id, fa_path)| {
        run_abundance_script(
            &qual_dir, sample_id, fa_path,
            &qual_db, &classify_path, &args.level, 0.0, &scripts_dir,
        );
    });

    // ══════════════════════════════════════════════════════════════════════════
    // PHASE 3 — Species identification filter
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 3: Species identification filter ════════════════════════");

    digest_outputs.par_iter().for_each(|(sample_id, _)| {
        let outdir   = format!("{qual_dir}/{sample_id}");
        let xls_path = format!("{outdir}/{sample_id}.xls");
        match &filter_mode {
            FilterMode::GScore(n) => run_gscore_filter(sample_id, &xls_path, &outdir, *n),
            FilterMode::ML(p)     => run_map2b_ml(sample_id, &xls_path, &outdir, *p),
        }
    });

    // ══════════════════════════════════════════════════════════════════════════
    // PHASE 4 — Build per-sample quantitative database
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 4: Build per-sample quantitative databases ══════════════");

    let mkdb_dir = format!("{output}/2.mkdb");
    fs::create_dir_all(&mkdb_dir).expect("Cannot create mkdb dir");

    let mkdb_script = find_sibling("scripts/CreatDB4AllLevel.py");
    let intersection_flag = if args.intersection { "--intersection" } else { "" };

    digest_outputs.par_iter().for_each(|(sample_id, fa_path)| {
        let smp_mkdb = format!("{mkdb_dir}/{sample_id}");
        fs::create_dir_all(&smp_mkdb)
            .unwrap_or_else(|e| panic!("Cannot create {smp_mkdb}: {e}"));

        let reads_list = format!("{smp_mkdb}/reads.list");
        {
            let mut f = BufWriter::new(
                File::create(&reads_list).expect("Cannot create reads.list"),
            );
            writeln!(f, "{sample_id}\t{fa_path}").unwrap();
        }

        let pred_result = format!("{qual_dir}/{sample_id}/pred.result");
        if !std::path::Path::new(&pred_result).exists() {
            eprintln!("WARN [{sample_id}]: pred.result not found — skipping mkdb");
            return;
        }

        let out_prefix = format!("{smp_mkdb}/{sample_id}");
        eprintln!("INFO [{sample_id}]: CreatDB4AllLevel.py");

        let mut cmd = std::process::Command::new("python3");
        cmd.args([
            mkdb_script.to_str().unwrap_or("scripts/CreatDB4AllLevel.py"),
            "-d", &database,
            "-c", &classify_path,
            "-p", &pred_result,
            "-o", &out_prefix,
            "-s", &args.cutoff.to_string(),
            "-n", "m",
            "-l", &args.level,
        ]);
        if !intersection_flag.is_empty() { cmd.arg(intersection_flag); }

        match cmd.status() {
            Ok(s) if s.success() =>
                eprintln!("INFO [{sample_id}]: mkdb done → {out_prefix}.marisa"),
            Ok(s) => eprintln!("WARN [{sample_id}]: CreatDB4AllLevel.py exited {s}"),
            Err(e) => eprintln!("WARN [{sample_id}]: Cannot run CreatDB4AllLevel.py: {e}"),
        }
    });

    // ══════════════════════════════════════════════════════════════════════════
    // PHASE 5 — Quantitative abundance
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 5: Quantitative abundance ═══════════════════════════════");

    let quan_dir     = format!("{output}/3.quan");
    let coverage_dir = format!("{output}/3.quan/.cov");
    fs::create_dir_all(&quan_dir)    .expect("Cannot create quan dir");
    fs::create_dir_all(&coverage_dir).expect("Cannot create coverage dir");

    let quan_script = scripts_dir.join("CalculateRelativeAbundance_Single2bEnzyme.py");

    digest_outputs.par_iter().for_each(|(sample_id, _)| {
        let mkdb_prefix = format!("{mkdb_dir}/{sample_id}/{sample_id}");
        let reads_list  = format!("{mkdb_dir}/{sample_id}/reads.list");
        let mkdb_marisa = format!("{mkdb_prefix}.marisa");

        if !std::path::Path::new(&mkdb_marisa).exists() {
            eprintln!("WARN [{sample_id}]: {mkdb_marisa} not found — skipping quan");
            return;
        }

        for (out_dir, ct) in [(&quan_dir, args.cov_thresh), (&coverage_dir, 0.0_f64)] {
            fs::create_dir_all(out_dir)
                .unwrap_or_else(|e| panic!("Cannot create {out_dir}: {e}"));
            eprintln!("INFO [{sample_id}]: quan ct={ct} → {out_dir}");
            match std::process::Command::new("python3")
                .args([
                    quan_script.to_str().unwrap_or("scripts/CalculateRelativeAbundance_Single2bEnzyme.py"),
                    "-d", &mkdb_prefix,
                    "-c", &classify_path,
                    "-l", &reads_list,
                    "-t", &args.level,
                    "-o", out_dir,
                    "-p", "1",
                    "-ct", &ct.to_string(),
                ])
                .status()
            {
                Ok(s) if s.success() => {},
                Ok(s) => eprintln!("WARN [{sample_id}]: quan (ct={ct}) exited {s}"),
                Err(e) => eprintln!("WARN [{sample_id}]: quan error: {e}"),
            }
        }
    });

    // ══════════════════════════════════════════════════════════════════════════
    // PHASE 6 — Merge results across all samples
    // ══════════════════════════════════════════════════════════════════════════
    eprintln!("\n═══ Phase 6: Merging results ══════════════════════════════════════");

    let stat_dir = format!("{output}/4.stat");
    fs::create_dir_all(&stat_dir).expect("Cannot create stat dir");

    let abd_list_path = format!("{stat_dir}/abd.list");
    {
        let mut f = BufWriter::new(File::create(&abd_list_path).expect("Cannot create abd.list"));
        for (sample_id, _) in &digest_outputs {
            writeln!(f, "{sample_id}\t{quan_dir}/{sample_id}/{sample_id}.xls").unwrap();
        }
    }

    let cov_list_path = format!("{stat_dir}/cov.list");
    {
        let mut f = BufWriter::new(File::create(&cov_list_path).expect("Cannot create cov.list"));
        for (sample_id, _) in &digest_outputs {
            writeln!(f, "{sample_id}\t{coverage_dir}/{sample_id}/{sample_id}.xls").unwrap();
        }
    }
    eprintln!("INFO: abd.list → {abd_list_path}");
    eprintln!("INFO: cov.list → {cov_list_path}");

    let run_merge = |script_name: &str, extra_args: &[&str]| {
        let script = find_sibling(&format!("scripts/{script_name}"));
        let mut cmd = std::process::Command::new("python3");
        cmd.arg(script.to_str().unwrap_or(script_name));
        cmd.args(extra_args);
        eprintln!("INFO: {script_name}");
        match cmd.status() {
            Ok(s) if s.success() => {},
            Ok(s)  => eprintln!("WARN: {script_name} exited {s}"),
            Err(e) => eprintln!("WARN: cannot run {script_name}: {e}"),
        }
    };

    run_merge("MergeProfilesFromMultipleSamples.py", &[
        "-l", &abd_list_path,
        "-o", &format!("{stat_dir}/Abundance.tsv"),
    ]);
    run_merge("MergeCoverageFromMultipleSamples.py", &[
        "-i", &cov_list_path,
        "-o", &format!("{stat_dir}/Coverage.tsv"),
    ]);

    let digest_list_path = format!("{stat_dir}/digest.list");
    {
        let mut f = BufWriter::new(File::create(&digest_list_path).expect("Cannot create digest.list"));
        for (sample_id, fa_path) in &digest_outputs {
            writeln!(f, "{sample_id}\t{fa_path}").unwrap();
        }
    }
    run_merge("dige_stat.py", &[
        "-i", &digest_dir,
        "-l", &digest_list_path,
        "-e", &args.enzyme,
        "-o", &format!("{stat_dir}/Tags.tsv"),
    ]);

    let metadata = format!("{db_parent}/metadata.tsv.gz");
    if std::path::Path::new(&metadata).exists() {
        run_merge("anno_abund_processor.py", &[
            "-i", &format!("{stat_dir}/Abundance.tsv"),
            "-d", &metadata,
            "-o", &stat_dir,
        ]);
    }

    eprintln!("\n═══ Complete ═══════════════════════════════════════════════════════");
    eprintln!("INFO: Digest       → {digest_dir}/");
    eprintln!("INFO: Final output → {stat_dir}/");
    eprintln!("INFO:   Abundance.tsv  — relative abundance (sum=1, species with coverage ≥ {})", args.cov_thresh);
    eprintln!("INFO:   Coverage.tsv   — coverage for all detected species");
    eprintln!("INFO:   Tags.tsv       — digest statistics");
}
