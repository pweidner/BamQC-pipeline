# BamQC-pipeline 🌻

Snakemake workflow for automated BAM alignment statistics intended to run on the BIH cluster using SLURM scheduler for parallelization.

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
snakemake --config ref=hg38 data_location=/path/to/bamdir output_location=/path/to/outputdir --profile workflow/profiles
```
When first executed, the pipeline automatically installs all tools from a conda qc_env.yaml, which can be reused for any other runs.

`ref` -> Currently specifies the path to my personal reference files I use for alignments, should be updated to group references soon.

Find a summary of the current output generated [here](workflow/Output.md).

## 4. 💂‍♂️ Authors 

- [Patrick Weidner](https://github.com/pweidner)
