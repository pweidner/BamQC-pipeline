# BamQC-pipeline 🚀

A fast, reproducible **Snakemake** workflow for **QC of BAMs** — alignment stats, binned-coverage metrics, library complexity, and optional **Ashley’s QC**.

**Highlights**
- 🧬 **Alfred** alignment + error metrics  
- 🧊 **Binned coverage** → entropy, spikiness, GC bias  
- 📈 **preseq** library complexity curves  
- 🧠 **Ashley’s QC**: auto-merge existing labels/features or predict  
- 🖼️ **Plots**: per-library PDFs + run summary  
- 🗂️ **Clean outputs** with consistent join key `Library`

---

## ⚙️ 1. Installation

```bash
git clone https://github.com/pweidner/BamQC-pipeline.git
cd BamQC-pipeline
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

---

## 🧠 5. Join keys

- Canonical key across all merges: **`Library`** (sample ID from BAM)
- `metadata/library_map.tsv` maps:
  ```
  cell                              Library
  A5455_L2_1001.sort.mdup.bam      OP-BB10-T_A5455_L2_1001
  ...
  ```
- `ashleys/prediction/prediction.norm.tsv` → `Library, ashleys_prediction, ashleys_probability, ashleys_sample`
- `ashleys/features.norm.tsv` → `Library, [feature columns…]`

---

## 🧬 6. Output Metrics

| **Category** | **Representative Columns** | **Description** |
|---------------|-----------------------------|-----------------|
| **Library Metadata** | `Library`, `Sample`, `sample_name`, `Base` | Unique identifiers for libraries/samples. Used as join key across all tables. |
| **Read Quality Filtering** | `#QCFail`, `QCFailFraction`, `#DuplicateMarked`, `DuplicateFraction`, `#Unmapped`, `UnmappedFraction` | Counts and fractions of reads filtered due to QC failure, duplication, or unmapping. Computed from **Alfred** metrics. |
| **Alignment Summary** | `#Mapped`, `MappedFraction`, `#MappedRead1`, `#MappedRead2`, `RatioMapped2vsMapped1`, `#MappedForward`, `MappedForwardFraction`, `#MappedReverse`, `MappedReverseFraction` | Alignment distribution across reads and orientations. |
| **Alignment Types** | `#SecondaryAlignments`, `SecondaryAlignmentFraction`, `#SupplementaryAlignments`, `SupplementaryAlignmentFraction`, `#SplicedAlignments`, `SplicedAlignmentFraction` | Multi-mapping, split, or spliced alignment metrics. |
| **Pairing & Proper Pairs** | `#Pairs`, `#MappedPairs`, `MappedPairsFraction`, `#MappedSameChr`, `MappedSameChrFraction`, `#MappedProperPair`, `MappedProperFraction` | Mate-pair statistics, intra-chromosomal mapping, and pairing quality. |
| **Reference / Coverage Stats** | `#ReferenceBp`, `#ReferenceNs`, `#AlignedBases`, `CoveredBp`, `FractionCovered`, `MedianCoverage`, `SDCoverage` | Total reference bases, aligned coverage, and coverage spread across genome. |
| **Base-level Alignment Accuracy** | `#MatchedBases`, `MatchRate`, `#MismatchedBases`, `MismatchRate`, `#DeletionsCigarD`, `DeletionRate`, `#InsertionsCigarI`, `InsertionRate`, `ErrorRate` | Alignment accuracy and indel rates, derived from **Alfred error metrics**. |
| **Clipping and Context** | `#SoftClippedBases`, `SoftClipRate`, `#HardClippedBases`, `HardClipRate`, `HomopolymerContextDel`, `HomopolymerContextIns` | Soft/hard clipping frequency and homopolymer indel contexts. |
| **Insert Size / Layout** | `DefaultLibraryLayout`, `MedianReadLength`, `MedianInsertSize_x`, `MedianInsertSize_y` | Median read and fragment lengths; layout (PE vs SE). |
| **Mapping Quality** | `MedianMAPQ_x`, `MedianMAPQ_y`, `p_mq` | Distribution and fraction of high-MAPQ reads. |
| **Binned Coverage Metrics** | `n.bins`, `avg.binsize`, `total.read.count`, `avg.read.count`, `spikiness`, `entropy`, `coverage_gini`, `coverage_cv`, `coverage_mad`, `coverage_sd` | **bedtools-binned** metrics characterizing coverage uniformity and noise. Entropy/spikiness quantify library evenness and strand-specific artifacts. |
| **GC Bias & Fold80** | `fold80_penalty`, `gc_pearson_r` | GC-coverage correlation and Fold-80 penalty (uniformity measure). |
| **Library Complexity (preseq)** | `preseq_distinct_at_observed`, `preseq_saturation` | Predicted complexity curve and saturation fraction at observed sequencing depth. |
| **Coverage Depth Fractions** | `pct_ge_1x`, `pct_ge_10x`, `pct_ge_30x` | Percentage of reference bases covered ≥ 1×, 10×, 30×. |
| **Ashley’s QC Predictions** | `prediction`, `probability` | Auto-classification of library quality via **Ashley’s QC** model or feature merge. |
| **Mosaicatcher Counts (good reads)** | `p_unmap`, `p_map`, `p_supp`, `p_dup`, `p_read2`, `p_good` | Fractions of reads in each category from `counts.info` — representing filtered “good” reads in downstream plots. |
| **Binned Window Weights** | `W10_5.0mb` → `W100_0.2mb`, `total_*` | Fractional read coverage across genome windows (multi-scale binning). Useful for CNV/spikiness consistency QC. |
````

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
