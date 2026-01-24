# BamQC-pipeline 🚀

A fast, reproducible **Snakemake** workflow for **QC of BAMs** — alignment stats, binned-coverage metrics, library complexity, and optional **Ashley’s QC**.

**Highlights**
- 🧬 **Alfred** alignment + error metrics  
- 🧊 **Bedtools** → binned coverage entropy, spikiness, gini, GC bias
- 📈 **Preseq** library complexity curves  
- 🧠 **ASHLEYS QC**: auto-merge existing labels/features or predict  
- 🖼️ **Plots**: per-library PDFs + run summary  
- 🗂️ **Clean outputs** with consistent join key `Library`

---

## ⚙️ 1. Installation

Clone the repo into your work folder:

```bash
git clone https://github.com/pweidner/bamqc-pipeline.git
cd bamqc-pipeline
```

You will need the same minimal snakemake env as for running mosaicatcher (you can also use that if you already have it, then just activate):

```bash
mamba create -n snakemake snakemake=7.32.0 -y
conda activate snakemake
```
> Tools used by rules are auto-installed from `workflow/envs/*.yaml` on first run.

---

## 🧾 2. Configuration (`config/config.yaml`)

```yaml
ref: hg38
reference_path: /ref/dir               # contains hg38.fa (+.fai)
data_location:  /path/to/input         # FLAT: *.sort.mdup.bam; HIER: <sample>/bam/*.sort.mdup.bam
output_location: /path/to/output
window: 200000
plot: true

bam_ext: ".sort.mdup.bam"
tmp_dir: /tmp

ashleys:
  enabled: true
  bin: /abs/path/ashleys-qc/bin/ashleys.py
  model_path: /abs/path/ashleys-qc/models/svc_default.pkl
  win_sizes: [5000000,2000000,1000000,800000,600000,400000,200000]
  threads: 32
  mem_mb: 200000
  conda_env: envs/ashleys.yaml
```

**Discovery modes**
- **FLAT**: `data_location/*.sort.mdup.bam`
- **HIER**: `data_location/<SAMPLE>/bam/*.sort.mdup.bam`  
  → auto-detected.

---

## ▶️ 3. Run

```bash
snakemake --config data_location=/data/runA output_location=/data/runA/bamqc --profile workflow/profiles --keep-going
```

---

## 📂 4. Output structure

```
output_location/
├── final_qc.tsv                         # 🧩 Alfred + counts-based + counts.info + preseq QC + Ashley’s columns
│
├── results/                             # 📊 Per tool deliverables
│   ├── final_qc.tsv
│   ├── preseq_metrics.tsv               # Preseq summary stats across libraries
│   └── alignment_summary_metrics.tsv    # Alfred summary across libraries
│
├── metadata/
│   └── library_map.tsv                  # cell <-> Library mapping for sanity checks
│
├── stats-by-lib/
│   └── {Library}.qc.tsv.gz              # per-lib Alfred output
│
├── binned/
│   └── {Library}.bins.tsv.gz            # windowed counts
│
├── qc-from-bins/
│   └── {Library}.counts_qc.tsv          # entropy/spikiness/GC-bias metrics
│
├── preseq/
│   └── {Library}.lc.tsv                 # library complexity curves
│
├── ashleys/
│   ├── features.tsv                     # merged/computed Ashley features
│   ├── features.norm.tsv                # features keyed by Library
│   └── prediction/
│       ├── prediction.tsv               # merged labels or predictions
│       └── prediction.norm.tsv          # normalized to Library
│
└── plots/
    ├── per-lib-qc/{Library}.qc.pdf      # optional per-lib PDF
    └── run_summary.pdf                  # run/cohort summary plots
```

## 🧬 Output metrics (what they mean + how to read them)

This pipeline produces per-library QC summaries in two main tables:

- **`final_qc.tsv`** — core QC metrics derived from **Alfred**, **bin-wise coverage**, and **preseq**
- **`final_qc_with_ashleys.tsv`** — `final_qc.tsv` plus **Ashley’s QC predictions** and **Ashley feature vectors**

All non-identifier columns are prefixed by their producing tool to make provenance explicit.

---

## 1️⃣ Identifiers (no prefix)

| Column | Description |
|------|-------------|
| `Library` | Unique library identifier used throughout the pipeline (e.g. `DRUG-CDTR-P1AZA-C_A5573_L1_i301.sort.mdup`). |
| `Sample` | Sample / condition identifier grouping multiple libraries (e.g. `DRUG-CDTR-P1AZA-C`). |

---

## 2️⃣ Alfred alignment & BAM QC (`alf_*`)

Metrics derived from **`alfred qc`**, summarizing mapping, alignment accuracy, and coverage statistics.

### Read filtering & mapping

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_qcfail_n` | QC-failed reads | High values indicate poor read quality. |
| `alf_qcfail_frac` | Fraction QC-failed | >0.05 often indicates a problematic library. |
| `alf_duplicate_marked_n` | Duplicate reads | |
| `alf_duplicate_frac` | Duplicate fraction | **High = low complexity or over-sequencing.** |
| `alf_unmapped_n` | Unmapped reads | |
| `alf_unmapped_frac` | Fraction unmapped | High values may indicate contamination or wrong reference. |
| `alf_mapped_n` | Mapped reads | |
| `alf_mapped_frac` | Fraction mapped | Healthy libraries are typically high (>0.8). |

### Read balance & orientation

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_mapped_read1_n`, `alf_mapped_read2_n` | Read1 / Read2 mapped counts | |
| `alf_mapped2_vs_mapped1_ratio` | Read2 / Read1 ratio | **≈1.0 expected** for paired-end data. |
| `alf_mapped_forward_frac` | Fraction forward strand | |
| `alf_mapped_reverse_frac` | Fraction reverse strand | **≈0.5 / 0.5 expected** unless protocol-biased. |

### Alignment types

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_secondary_alignments_frac` | Secondary alignments | Elevated values indicate multi-mapping / repeats. |
| `alf_supplementary_alignments_frac` | Supplementary alignments | Can indicate SVs, chimeras, or mapping artifacts. |
| `alf_spliced_alignments_frac` | Spliced alignments | Typically low for DNA; higher for RNA or odd mapping. |

### Pairing & concordance

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_mapped_pairs_frac` | Fraction of mapped read pairs | |
| `alf_mapped_same_chr_frac` | Pairs on same chromosome | Low values may indicate discordant mapping. |
| `alf_mapped_proper_pair_frac` | Properly paired reads | **High = good library structure.** |

### Alignment accuracy & errors

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_match_frac` | Matched bases / aligned bases | **Closer to 1 = higher accuracy.** |
| `alf_mismatch_frac` | Mismatched base fraction | |
| `alf_deletion_frac`, `alf_insertion_frac` | Indel rates | Elevated values can indicate mapping or chemistry artifacts. |
| `alf_error_frac` | Aggregate alignment error rate | Lower is better. |

### Clipping & context

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_soft_clip_frac` | Soft-clipped bases | High = adapters, short inserts, or mapping difficulty. |
| `alf_hard_clip_frac` | Hard-clipped bases | Usually near zero. |
| `alf_homopolymer_context_del/ins` | Indels in homopolymers | Elevated values indicate systematic indel artifacts. |

### Read length, coverage & MAPQ

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `alf_read_length_med` | Median read length | |
| `alf_insert_size_med` | Median insert size | |
| `alf_mapq_med` | Median mapping quality | **Higher = more confident mapping.** |
| `alf_coverage_med` | Median coverage | |
| `alf_coverage_sd` | Coverage standard deviation | High = uneven coverage. |
| `alf_covered_frac` | Fraction of reference covered | Low = sparse library. |

---

## 3️⃣ Bin-wise coverage metrics (`bin_*`)

Computed from **fixed-size genome windows** using `bedtools coverage -counts` and summarized in `qc_from_counts.py`.

### Basic bin descriptors

| Column | Meaning |
|------|--------|
| `bin_n_bins` | Number of windows used. |
| `bin_avg_binsize` | Mean window size (bp). |
| `bin_total_read_count` | Total reads across all bins. |
| `bin_avg_read_count` | Mean reads per bin. |

### Coverage uniformity & signal shape

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `bin_entropy` | Shannon entropy of bin counts | **Higher = more even coverage.** |
| `bin_spikiness` | Local coverage jaggedness | **Higher = noisier / uneven signal.** |
| `bin_gini` | Gini index of coverage | **0 = uniform, higher = uneven.** |
| `bin_cv` | Coefficient of variation | |
| `bin_mad` | Median absolute deviation | Robust variability measure. |
| `bin_sd` | Standard deviation | |

### Uniformity & GC bias

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `bin_fold80` | Fold-80 penalty | **≈1 ideal**, higher = worse uniformity. |
| `bin_gc_r` | GC–coverage correlation | Large magnitude = GC bias. |

### Depth thresholds

| Column | Meaning |
|------|--------|
| `bin_pct_ge_1x` | Fraction of bins with ≥1 read. |
| `bin_pct_ge_10x` | Fraction with ≥10 reads. |
| `bin_pct_ge_30x` | Fraction with ≥30 reads. |

---

## 4️⃣ Library complexity (preseq) (`preseq_*`)

Derived from **`preseq lc_extrap`**, estimating how many *unique* DNA fragments are present.

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `preseq_distinct_at_observed` | Expected number of distinct fragments at observed depth | Higher = more complex library. |
| `preseq_saturation` | Distinct / total reads at observed depth | **0 = highly duplicated**, **1 = highly complex**. |

**Rule of thumb**
- `preseq_saturation ≈ 1` → sequencing deeper will still yield new information  
- `preseq_saturation ≈ 0` → sequencing deeper mostly yields duplicates  

---

## 5️⃣ Ashley QC predictions (`ash_*`)

Generated by **ashleys-qc**, integrating coverage patterns and strand balance.

| Column | Meaning | Interpretation |
|------|--------|----------------|
| `ash_label` | Predicted QC class (model-specific) | |
| `ash_prob` | Prediction confidence | Values near 0.5 = ambiguous. |
| `ash_cell` | BAM identifier used by Ashley | |
| `ash_sample` | Sample label used by Ashley | |

---

## 6️⃣ Ashley feature vectors (`ash_*`)

Multi-scale **Watson-strand bin features** and read category fractions.

### Window-bin distributions

For each window size (`5mb`, `2mb`, `1mb`, `0_8mb`, `0_6mb`, `0_4mb`, `0_2mb`):

| Column pattern | Meaning |
|--------------|--------|
| `ash_w10_*` … `ash_w100_*` | Fraction of windows falling into Watson% bins (0–10%, …, 90–100%). |
| `ash_total_*` | Total fraction across bins (≈1.0 if normalized). |

**Interpretation**
- Smooth, balanced distributions indicate stable coverage.
- Skewed distributions can indicate CNVs, strand imbalance, or technical artifacts.

### Mapping & “good read” fractions

| Column | Meaning |
|------|--------|
| `ash_p_unmap` | Fraction unmapped. |
| `ash_p_map` | Fraction mapped. |
| `ash_p_supp` | Fraction supplementary. |
| `ash_p_dup` | Fraction duplicate. |
| `ash_p_mq` | Fraction passing MAPQ filter. |
| `ash_p_read2` | Fraction read2. |
| `ash_p_good` | Fraction of reads passing all Ashley filters (usable signal). |

---

## 🔍 Practical interpretation summary

- **Low complexity** → high `alf_duplicate_frac`, low `preseq_saturation`
- **Uneven coverage** → high `bin_spikiness`, `bin_gini`, `bin_fold80`
- **GC bias** → large `|bin_gc_r|`
- **Mapping problems** → low `alf_mapped_frac`, low `alf_mapq_med`
- **Ashley disagreement** → `ash_prob low` → inspect manually

---



## Notes
- **Entropy** and **spikiness** reflect coverage evenness (low entropy or high spikiness = uneven).
- **Fold80 penalty** follows the Picard metric (ideal = 1, higher = less uniform).
- **Preseq** metrics allow extrapolation of unique reads vs sequencing depth.
- **Ashley’s QC** integrates pretrained classification of Strand-seq libraries by coverage pattern and W→C balance.
- **Mosaicatcher** fractions (`p_good`, etc.) summarize the final usable subset for downstream analyses like count plots and phasing.

---

## 📚 8. Citations

- **Alfred** — Rausch *et al.*, *Genome Res* (2019)  
- **preseq** — Daley & Smith, *Bioinformatics* (2013)
- **ASHLEYS** - Gros *et al.*, *Bioinformatics* (2021)
- **bedtools** — Quinlan & Hall, *Bioinformatics* (2010)  
- **Snakemake** — Köster & Rahmann, *Bioinformatics* (2012)

---

## 💡 9. Roadmap

- contamination checks  
- More GC/coverage plots (Lorenz, violin)  
- Optional HTML dashboard  

---

Happy QC’ing! 🧪✨
