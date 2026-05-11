# VIP2B
VIrome Profiler with type IIB restriction sites for WMS data

---

## Installation
```bash
git clone https://github.com/YiyanYang0728/VIP2B.git
cd VIP2B
# Use the provided conda environment
conda create -n VIP2B -c conda-forge --file requirement.txt
conda activate VIP2B
# If you have mamba, run
mamba create -n VIP2B -c conda-forge --file requirement.txt
mamba activate VIP2B
```
Make sure the conda environment of VIP2B has been activated by running the above command before you run VIP2B everytime.

Build an executable software:
```bash
cargo build --release
# Add the executive software to PATH:
export PATH=$PATH:$PWD/target/release
# Or copy it to your bin directory: cp target/release/VIP2B <your bin path>
```

Download databases:
```bash
mkdir -p database
# this may take a while
python scripts/DownloadDB.py -l config/def_db.list -d database
```

---

## Quick start & test
```bash
# We prepared a pair-end and a single-end sequencing test data:
cd example
mkdir -p data/
wget -t 3 -O data/test_seq.R1.fq.gz https://zenodo.org/records/19700476/files/test_seq.R1.fq.gz
wget -t 3 -O data/test_seq.R2.fq.gz https://zenodo.org/records/19700476/files/test_seq.R2.fq.gz
wget -t 3 -O data/test_seq.fq.gz https://zenodo.org/records/19700476/files/test_seq.fq.gz
db=$PWD/../database # or <your full path for the database>
VIP2B -i data.list -p 4 -d ${db}/8Enzyme -o test_result
```

---

## Sample list format (`-i`)

A tab-separated file — one sample per line. Lines beginning with `#` are ignored.
Single-end and paired-end samples can be mixed.
```
# single-end (2 columns)
sample1    /data/sample1.fastq.gz
# paired-end (3 columns)
sample2    /data/sample2_R1.fastq.gz    /data/sample2_R2.fastq.gz
```

---

## parameters
```
VIP2B -h
Unified 2bRAD virome profiler (digest + abundance)

Usage: VIP2B [OPTIONS] -i <FILE>

Options:
  -i <FILE>             Sample list TSV (2-col SE or 3-col PE; # lines ignored) Format: sample_id <TAB> reads.fastq.gz [<TAB> reads_R2.fastq.gz]
  -o <DIR>              Output directory [default: VIP2B_result in current directory] [default: VIP2B_result]
  -l <LEVEL>            Taxonomy level [Class|Order|Family|Genus|Species] [default: Species] [possible values: Class, Order, Family, Genus, Species]
  -e <ENZYME[,...]>     Enzyme(s): comma-separated names or "all". Available: AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I The default 8-enzyme set matches the standard database [default: AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I]
  -d <PREFIX>           Database prefix (without extension, e.g. /path/to/database/8Enzyme) [default: database/8Enzyme relative to the vip2b binary] [default: database/8Enzyme]
  -p <N>                Number of parallel processes/threads (more may require more memory) [default: 1]
  -t <G5|M0.5>          Threshold for species identification. G
                        : G-score filter applied in Rust (e.g. G2, G5). M{p}: MAP2B_ML.py is called automatically after the abundance step [default: M0.5]
  -c <N>                Cut-off for database building (mkdb step only; not applied in Rust) [default: 30000]
  -f <FLOAT>            Threshold for coverage filtering: sequenced_tags / theoretical_tags [default: 0.6]
      --intersection    Use intersection of tags between genomes (default: union). Applies to mkdb step only; stored here for pipeline compatibility
      --batch-size <N>  Reads per parallel batch in digest phase [default: 100000]
  -h, --help            Print help (see more with '--help')
  -V, --version         Print version
```

Parameter explanation:

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Sample list TSV (2-col SE or 3-col PE; required) | — |
| `-o` | Output directory | `./VIP2B_result` |
| `-l` | Taxonomy level: `Class` / `Order` / `Family` / `Genus` / `Species` | `Species` |
| `-e` | Enzyme(s), comma-separated, or `all` for all 8 enzymes | `AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I` |
| `-d` | Database prefix — path without extension (e.g. `/db/8Enzyme`) | `../database/8Enzyme` |
| `-p` | Number of parallel processes/threads | `1` |
| `-t` | Species ID threshold: `G{n}` (G-score) or `M{p}` (ML probability) | `M0.5` |
| `-c` | Cut-off for per-sample DB building (mkdb step) | `30000` |
| `-f` | Coverage threshold for abundance output (sequenced / theoretical tags) | `0.6` |
| `--intersection` | Use tag intersection instead of union in mkdb step | off |
| `--batch-size` | Reads per Rust parallel batch (tuning only) | `100000` |

Available enzymes (16 total):  
`AlfI AloI BaeI BcgI BplI BsaXI BslFI Bsp24I CjeI CjePI CspCI FalI HaeIV Hin4I PpiI PsrI`

### Species ID threshold (`-t`)

| Value | Behaviour |
|-------|-----------|
| `G2`, `G5`, … | G-score > n — `gscore_filter.py` called automatically |
| `M0.1`, `M0.5`, … | ML probability > p — `MAP2B_ML.py` called automatically (requires XGBoost) |

---

## Pipeline phases and output structure

The pipeline produces **5 output folders**, matching the original VIP2B exactly:

```
{output}/
├── 0.dige/                         Phase 1 — Tag extraction (Rust)
│   ├─ {sample}.fa.gz              Extracted 2bRAD tags (gzip FASTA)
│   └── {sample}/
│       └── {sample}.dige.stat.xls  Digest statistics (input reads / extracted tags / rate)
│
├── 1.qual/                         Phase 2+3 — Qualitative abundance + species filter
│   └── {sample}/
│       ├── {sample}.xls            Per-sample abundance (all candidate species, ct=0)
│       ├── {sample}.GCF_detected.xls  GCF-level coverage
│       └── pred.result             Species that passed the -t filter (input to mkdb)
│
├── 2.mkdb/                         Phase 4 — Per-sample quantitative database
│   └── {sample}/
│       ├── {sample}.marisa         Per-sample tag trie (only confirmed species)
│       ├── {sample}.stat.xls       Theoretical tag counts for confirmed species
│       └── reads.list              Sample→FASTQ mapping for this sample
│
├── 3.quan/                         Phase 5 — Quantitative abundance
│   └── {sample}/
│       ├── {sample}.xls            Abundance filtered by -f threshold (→ Abundance.tsv)
│       └── {sample}.GCF_detected.xls
│
└── 4.stat/                         Merged results (final output)
    ├── Abundance.tsv               Relative abundance per species
    ├── Coverage.tsv                Taxonomic coverage per species
    ├── ...                         Other inferred profiles. Please see details in `Final output files`
    └── Tags.tsv                    Digest statistics across all samples
```

### Final output files

| File | Description |
|------|-------------|
| `4.stat/Abundance.tsv` | Relative abundance table — species × samples, values sum to 1 per sample. Only species with coverage ≥ `-f` threshold. |
| `4.stat/Coverage.tsv` | Taxonomic coverage table — species × samples, value = sequenced_tags / theoretical_tags. All detected species regardless of coverage. |
| `4.stat/Phenotype.tsv` | Viral phenotype profiles — jumbo phage and lytic phage percentage per sample. |
| `4.stat/viral_function/Uniref90.tsv` | Gene family abundance table (Uniref90 annotations). |
| `4.stat/viral_function/cluster.tsv` | Gene cluster abundance table. |
| `4.stat/host_taxonomy/<taxonomy>_abund.tsv` | Virus host abundance at a certain taxonomy level. |

> **Note:** `4.stat/Phenotype.tsv`, `4.stat/viral_function/`, and `4.stat/host_taxonomy/` are produced by `anno_abund_processor.py` and require `metadata.tsv.gz` to be present in the database directory. 

---

## File structure

```
VIP2B/
├── Cargo.toml                  Rust workspace definition
├── README.md
├── requirement.txt             Python dependency list
├── config/
│   ├── def_db.list             Database file list and Zenodo URLs for auto-download
│   ├── requirement.txt         Python dependencies (conda environment spec)
│   └── XGB_none_0238.pkl       XGBoost model for -t M species ID filtering
├── example/
│   └── data.list               Example sample list
├── scripts/                    Python helper scripts (called automatically by VIP2B)
│   ├── CalculateRelativeAbundance_Single2bEnzyme.py   Phase 2 — qualitative abundance
│   ├── CreatDB4AllLevel.py         Build per-sample quantitative database
│   ├── quan_from_tags.py           marisa lookup from Rust-piped tag counts
│   ├── gscore_filter.py            G-score species ID filter (-t G mode)
│   ├── MAP2B_ML.py                 ML species ID filter (-t M mode)
│   ├── MergeProfilesFromMultipleSamples.py   produce Abundance.tsv
│   ├── MergeCoverageFromMultipleSamples.py   produce Coverage.tsv
│   ├── dige_stat.py                produce Tags.tsv
│   ├── anno_abund_processor.py     produce Phenotype.tsv, viral_function/, host_taxonomy/ (requires metadata.tsv.gz)
│   ├── DownloadDB.py               Auto-download default database from Zenodo
│   ├── convert_marisa.py           Utility: convert marisa → TSV (not used by pipeline)
│   ├── assess.py
│   ├── host_filter.py
│   ├── host_marisa_trie.build.py
│   ├── marisa_trie.build.py
│   └── sequence_digestion.py       Original Python digest (superseded by Rust)
└── vip2b/
    ├── Cargo.toml              Rust package config (builds target/release/VIP2B)
    └── src/
        └── main.rs             Full pipeline source
```

