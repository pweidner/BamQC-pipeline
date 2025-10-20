# BamQC-pipeline 🌻

A Snakemake workflow for generating **comprehensive QC metrics** from BAM files, combining:
- **[Alfred](https://github.com/tobiasrausch/alfred)** alignment and error statistics
- **Bin-based QC metrics** (entropy, spikiness, coverage uniformity, GC bias)
- **[preseq](http://smithlabresearch.org/software/preseq/)** library complexity estimates
- Automated plotting of per-sample QC reports

---


## 1. 📕 Installation

Clone the repository into your personal work (e.g. /fast/users/$USER/work/)

```bash
git clone https://github.com/pweidner/BamQC-pipeline.git
```
The pipeline currently looks like this:

```
BamQC-pipeline/            # Run the pipeline from here
│
├── config/
│   └── config.yaml        # Config parameters for user input
└── workflow/
    ├── envs/
    │   └── qc_env.yaml    # Conda environment housing all tools (automated installation upon first execution)
    ├── errors/            # Error files by rules (generated upon first execution)
    ├── logs/              # Log files by rules (generated upon first execution)
    ├── profiles/
    │   └── config.yaml    # BIH specific SLURM config
    ├── Snakefile
    └── utils.py    
    
```

To run the pipeline, you need to set up a conda environment to launch snakemake from

```bash
mamaba create -n snakemake snakemake=7.32.0
```

## 2. 🛑 Run the pipeline

To launch the pipeline, simply cd into the folder containing the pipeline and activate the snakemake env:

```bash
cd /path/to/BamQC-pipeline
conda activate snakemake
```
Once inside the repo, you can launch any number of runs from here on I/O folders (e.g. from within a dedicated screen session):

```bash
snakemake --config ref=/path/to/hg38 data_location=/path/to/bamdir output_location=/path/to/outputdir --profile workflow/profiles
```
When first executed, the pipeline automatically installs all tools from a conda qc_env.yaml, which can be reused for any other runs.

`ref` -> Currently specifies the path to my personal reference files I use for alignments, should be updated to group references soon (without the .fa ending).


---

## 📂 Outputs

All outputs are written under the configured `output_location` (default: `results/`).

### 1. **Alfred QC**
- Path: `stats-by-lib/{sample}.qc.tsv.gz`
- Per-sample raw QC report from **alfred qc**.
- Includes:
  - Mapping stats (#Mapped, #Unmapped, fractions)
  - Duplicate rate
  - Error/mismatch rates
  - Insert size distribution
  - Base qualities
  - Coverage summary, etc.

### 2. **Alignment Summary**
- Path: `alignment_summary_metrics.tsv`
- Aggregated summary table across all samples (parsed from Alfred).
- One row per library.

### 3. **Binned Coverage**
- Path: `binned/{sample}.bins.tsv.gz`
- Read counts per genomic window (window size set via `config.yaml`).
- Generated with `samtools` + `bedtools coverage`.
- Used for downstream QC metrics.

### 4. **Counts-Based QC**
- Path: `qc-from-bins/{sample}.counts_qc.tsv`
- Per-sample QC metrics derived from windowed coverage:
  - **Entropy**: Shannon entropy of coverage distribution (higher = more uniform).
  - **Spikiness**: Bin-to-bin variability in coverage.
  - **Coverage CV / SD / MAD**: Variability measures.
  - **Coverage Gini**: Inequality of coverage distribution (0=perfectly uniform, 1=highly uneven).
  - **Fold-80 penalty**: Depth needed to bring 80% of bases to mean coverage.
  - **PCT≥X**: Fraction of bases covered at ≥1x, 10x, 30x (configurable).
  - **GC Pearson r**: Correlation between GC fraction and normalized coverage (GC bias).
  - Optional: insert size / MAPQ median (from Alfred), if available.

### 5. **Library Complexity (preseq)**
- Path: `preseq/{sample}.lc.tsv`
- Predicted distinct read counts at increasing sequencing depths.
- Used to assess library complexity / duplication.
- Reported as:
  - **preseq_distinct_at_observed**: Distinct reads at current depth.
  - **preseq_saturation**: Fraction of potential library complexity reached.

### 6. **Final Aggregated QC**
- Path: `final_qc.tsv`
- Merged table of **Alfred + counts-based + preseq metrics**.
- Master file for downstream analysis and plotting.

### 7. **Plots (optional, if `plot: true` in config)**
- Path: `plots/{sample}.qc.pdf`
- Multi-page PDF per sample containing:
  - **Summary page**: Key numbers (mapped fraction, duplicates, error rate, coverage, entropy, spikiness).
  - **Bar plots**: Key mapping and error fractions.
  - **Counts-based QC**: Entropy, spikiness, GC bias metrics.
  - **Preseq curve**: Predicted distinct vs. total reads (library complexity).

---

## ⚙️ Configuration (`config.yaml`)

Key settings:
- `data_location`: path to BAM files (`*.sort.mdup.bam`)
- `output_location`: results folder
- `ref`: reference genome name
- `reference_path`: path to reference FASTA
- `window`: bin size for coverage counts (e.g. 200000)
- `plot`: `true`/`false` — enable per-sample PDF reports

---

## 🗂 Example File Tree

Key settings:
- `data_location`: path to BAM files (`*.sort.mdup.bam`)
- `output_location`: results folder
- `ref`: reference genome name
- `reference_path`: path to reference FASTA
- `window`: bin size for coverage counts (e.g. 200000)
- `plot`: `true`/`false` — enable per-sample PDF reports

---

## 🗂 Example File Tree

results/
├── alignment_summary_metrics.tsv
├── final_qc.tsv
├── stats-by-lib/
│ └── sampleA.qc.tsv.gz
├── binned/
│ └── sampleA.bins.tsv.gz
├── qc-from-bins/
│ └── sampleA.counts_qc.tsv
├── preseq/
│ └── sampleA.lc.tsv
└── plots/
└── sampleA.qc.pdf

yaml
Code kopieren

---

## 📖 Metric Definitions (Quick Reference)

- **MappedFraction**: proportion of reads successfully mapped.
- **DuplicateFraction**: fraction of reads marked as PCR/optical duplicates.
- **ErrorRate**: base error rate across aligned reads.
- **Entropy**: diversity/uniformity of coverage across bins.
- **Spikiness**: relative variability in adjacent bin coverage.
- **Coverage CV / SD / MAD**: different measures of coverage spread.
- **Coverage Gini**: inequality index for coverage distribution.
- **Fold-80 Penalty**: ratio of mean coverage to coverage at 20th percentile of bins (lower = more uniform).
- **PCT≥X**: % of genome covered at or above X depth.
- **GC Pearson r**: correlation between GC content and coverage.
- **Preseq Saturation**: predicted fraction of unique molecules sequenced at current depth.

---

## 📑 Data Dictionary — `final_qc.tsv`

Each row is a library/sample. Columns come from:
- **Alfred** (`stats-by-lib/*.qc.tsv.gz` → parsed “ME” lines)
- **Counts-based QC** (`qc-from-bins/*.counts_qc.tsv`)
- **preseq** (`preseq/*.lc.tsv`)
- Merge occurs on **Library**

### Identifiers
| Column | Definition |
|---|---|
| `Library` | Library ID (derived from BAM filename prefix). |
| `Sample`  | Biological sample identifier parsed by Alfred. |

### Basic QC (Alfred)
| Column | Definition |
|---|---|
| `#QCFail` | Number of reads failing vendor/Illumina QC. |
| `QCFailFraction` | `#QCFail / total_reads`. |
| `#DuplicateMarked` | Reads flagged as duplicate (0x400). |
| `DuplicateFraction` | `#DuplicateMarked / total_reads`. |
| `#Unmapped` | Reads unmapped (0x4). |
| `UnmappedFraction` | `#Unmapped / total_reads`. |
| `#Mapped` | Reads mapped (primary). |
| `MappedFraction` | `#Mapped / total_reads`. |
| `#MappedRead1` | R1 mapped. |
| `#MappedRead2` | R2 mapped. |
| `RatioMapped2vsMapped1` | `#MappedRead2 / #MappedRead1`. |
| `#MappedForward` | Mapped on forward strand. |
| `MappedForwardFraction` | `#MappedForward / #Mapped`. |
| `#MappedReverse` | Mapped on reverse strand. |
| `MappedReverseFraction` | `#MappedReverse / #Mapped`. |
| `#SecondaryAlignments` | Secondary alignments (0x100). |
| `SecondaryAlignmentFraction` | `/ total_reads`. |
| `#SupplementaryAlignments` | Supplementary alignments (0x800). |
| `SupplementaryAlignmentFraction` | `/ total_reads`. |
| `#SplicedAlignments` | Spliced alignments (primarily RNA-seq; intron bridging). |
| `SplicedAlignmentFraction` | `/ #Mapped`. |

### Pairing / Insert stats (Alfred)
| Column | Definition |
|---|---|
| `#Pairs` | Total read pairs. |
| `#MappedPairs` | Pairs with both mates mapped. |
| `MappedPairsFraction` | `#MappedPairs / #Pairs`. |
| `#MappedSameChr` | Both mates on same chromosome. |
| `MappedSameChrFraction` | `/ #MappedPairs`. |
| `#MappedProperPair` | Properly paired by aligner (0x2). |
| `MappedProperFraction` | `/ #MappedPairs`. |
| `DefaultLibraryLayout` | SE/PE. |
| `MedianInsertSize` | Median template length (bp) for proper pairs. |

### Alignment & error profile (Alfred)
| Column | Definition |
|---|---|
| `#ReferenceBp` | Reference size considered (bp). |
| `#ReferenceNs` | ‘N’ bases in reference. |
| `#AlignedBases` | Total aligned query bases. |
| `#MatchedBases` | Matched aligned bases. |
| `MatchRate` | `#MatchedBases / #AlignedBases`. |
| `#MismatchedBases` | Mismatched aligned bases. |
| `MismatchRate` | `#MismatchedBases / #AlignedBases`. |
| `#DeletionsCigarD` | Count of deleted bases (CIGAR D). |
| `DeletionRate` | `/ #AlignedBases`. |
| `HomopolymerContextDel` | Fraction of deletions within homopolymers. |
| `#InsertionsCigarI` | Count of inserted bases (CIGAR I). |
| `InsertionRate` | `/ #AlignedBases`. |
| `HomopolymerContextIns` | Fraction of insertions within homopolymers. |
| `#SoftClippedBases` | Soft-clipped bases (CIGAR S). |
| `SoftClipRate` | `/ #AlignedBases`. |
| `#HardClippedBases` | Hard-clipped bases (CIGAR H). |
| `HardClipRate` | `/ #AlignedBases`. |
| `ErrorRate` | Typically `MismatchRate + InsertionRate + DeletionRate`. |
| `MedianReadLength` | Median read length (bp). |
| `MedianCoverage` | Alfred’s median per-base coverage. |
| `SDCoverage` | Standard deviation of per-base coverage. |
| `CoveredBp` | Number of reference bp with coverage ≥1×. |
| `FractionCovered` | `CoveredBp / (#ReferenceBp - #ReferenceNs)`. |
| `BpCov1ToCovNRatio` | Ratio of bp at coverage ≥1X to ≥N (N from Alfred). |
| `BpCov1ToCov2Ratio` | Ratio of bp ≥1X to ≥2X. |
| `MedianMAPQ` | Median mapping quality. |

### Binning metadata (counts input)
| Column | Definition |
|---|---|
| `n.bins` | Number of windows (bins) counted. |
| `avg.binsize` | Mean bin size (bp). *(== config `window` if constant)* |
| `total.read.count` | Sum of counts across all bins (filtered: no dups/secondary/supplementary; MAPQ≥10; chr1–22,X,Y). |
| `avg.read.count` | `total.read.count / n.bins`. |

### Counts-based QC (from `qc_from_counts.py`)
| Column | Definition |
|---|---|
| `entropy` | Shannon entropy of normalized bin coverage. Higher → more uniform. |
| `spikiness` | Adjacent-bin variability (e.g., mean absolute diff / mean). Higher → more spiky coverage. |
| `coverage_sd` | Standard deviation of (normalized) bin counts. |
| `coverage_cv` | Coefficient of variation = `coverage_sd / mean`. |
| `coverage_mad` | Median absolute deviation of (normalized) bin counts. |
| `coverage_gini` | Gini index of coverage inequality (0 uniform → 1 highly uneven). |
| `fold80_penalty` | `mean_coverage / P20_coverage` (depth ratio to bring 80% of genome to mean). |
| `pct_ge_1x` | % of bins with depth ≥1× (using normalized per-bin depth). |
| `pct_ge_10x` | % bins with depth ≥10×. |
| `pct_ge_30x` | % bins with depth ≥30×. |
| `gc_pearson_r` | Pearson correlation between per-bin GC fraction and normalized coverage (GC bias; negative → dropout at high GC). |
| `gc_slope` | Slope from linear fit: coverage ~ GC (optional; present if script outputs). |
| `gc_intercept` | Intercept of the GC linear fit (optional). |

> **Notes:**  
> • “Normalized coverage” in the table above refers to per-bin counts scaled (e.g., counts / mean counts) after optional bias trimming; check `qc_from_counts.py` for the exact normalization implemented.  
> • If `gc_slope`/`gc_intercept` are absent, only `gc_pearson_r` is reported.

### Library Complexity (preseq)
| Column | Definition |
|---|---|
| `preseq_distinct_at_observed` | Distinct reads (unique molecules) predicted at the observed sequencing depth (interpolated from `preseq` curve). |
| `preseq_total_reads_observed` | Observed total read count used for the above prediction. |
| `preseq_saturation` | `preseq_distinct_at_observed / preseq_total_reads_observed` (fraction of reads that are unique at current depth). |
| `preseq_curve_path` | Path to full preseq output table (`preseq/{sample}.lc.tsv`) for the extrapolation curve. |

---
## 📊 Citation

If you use this workflow, please cite:
- **Alfred**: Rausch T et al., *Genome Research* (2019).  
- **preseq**: Daley & Smith, *Bioinformatics* (2013).  
- **bedtools**: Quinlan & Hall, *Bioinformatics* (2010).  
- **Snakemake**: Köster & Rahmann, *Bioinformatics* (2012).

---

## 🚀 Future Extensions

Planned features:
- VerifyBamID / contamination checks
- Additional plots (coverage Lorenz curve, GC bias violin plots)
- Interactive dashboard (e.g. MultiQC-like HTML)

---

