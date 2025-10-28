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
snakemake   --config ref=hg38 reference_path=/ref            data_location=/data/runA            output_location=/data/runA/bamqc   --profile workflow/profiles   --rerun-incomplete --keep-going
```

---

## 📂 4. Output structure

All paths relative to your `output_location`.

```
output_location/
├── alignment_summary_metrics.tsv        # 🧬 Alfred summary across libraries
├── final_qc.tsv                         # 🧩 Alfred + counts-based + preseq QC
├── final_qc_with_ashleys.tsv            # 🧠 final_qc + Ashley’s columns
│
├── results/                             # 📊 deliverables
│   ├── final_qc.tsv
│   ├── final_qc_with_ashleys.tsv
│   └── alignment_summary_metrics.tsv
│
├── metadata/
│   └── library_map.tsv                  # cell <-> Library mapping
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
│   └── {Library}.lc.tsv                 # library complexity curve
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
    └── run_summary.pdf                  # optional run summary plot
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

## 🧬 6. Ashley’s integration

- If `*/cell_selection/labels.tsv` or similar exist → merged into `ashleys/prediction/prediction.tsv`
- If `*/predictions/ashleys_features.tsv` exist → merged into `ashleys/features.tsv`
- Otherwise, features + predictions are **computed** automatically  
- Everything is **normalized** to `Library` for merging

---

## 🖨️ 7. Quick examples

**Run full QC**
```bash
snakemake --config ref=hg38 data_location=/bams output_location=/out --profile workflow/profiles
```

**Generate only plots**
```bash
snakemake plots/run_summary.pdf --profile workflow/profiles
```

**Recompute Ashley features only**
```bash
snakemake ashleys/features.tsv --cores 16
```

---

## 📚 8. Citations

- **Alfred** — Rausch *et al.*, *Genome Res* (2019)  
- **preseq** — Daley & Smith, *Bioinformatics* (2013)
- **ASHLEYS** - Gros *et al.*, *Bioinformatics* (2021)
- **bedtools** — Quinlan & Hall, *Bioinformatics* (2010)  
- **Snakemake** — Köster & Rahmann, *Bioinformatics* (2012)

---

## 💡 9. Roadmap

- VerifyBamID / contamination checks  
- More GC/coverage plots (Lorenz, violin)  
- Optional HTML dashboard  

---

### ❤️ Contributing

Please open issues or PRs with:
- your `config.yaml`
- the `snakemake` command used
- relevant logs under `logs/` or `errors/`

Happy QC’ing! 🧪✨
