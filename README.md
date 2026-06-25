[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19854065.svg)](https://doi.org/10.5281/zenodo.19854065)

# Recurrent pseudoexon activation in COL4A5 intron 6 defines a therapeutically actionable hotspot

**Authors:** Hassan Saei, Béatrice Ardin, Nicolas Kaiser, Mouad Ouahmane, Olivier Gribouval, Vincent Moriniere, Florian Wopperer, Korbinian Riedhammer, Corinne Antignac, Michael Wiesener, Guillaume Dorval

This repository contains scripts and notebooks used to generate figures for the manuscript on deep intronic variation hotspot(s) in `COL4A5` responsible for X-linked Alport syndrome. 

The project is organized for reproducibility and clear separation between:
- AlphaGenome in-silico analyses (variant scoring, ISM saturation mutagenesis)
- Full-gene COL4A5 intron extraction and per-intron ISM
- Downstream postprocessing, motif enrichment, and summary plotting
- Complementary in-R motif GSEA cross-check
- Targeted RNA-seq sashimi visualization of the intron 6 hotspot

## Data availability
Targeted RNA-seq BAM files (+ `.bai` indexes) are archived on Zenodo:

**DOI:** [10.5281/zenodo.19854065](https://doi.org/10.5281/zenodo.19854065)

Download the archive and place the files in `targeted_rnaseq/bams/` (preserving filenames) before reproducing the sashimi figure. See `targeted_rnaseq/README.md` for sample design and ggsashimi commands.

## Repository Structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── coda/                              # Conda environment specs
│   ├── README.md
│   ├── alphagenome-env.yml            # Python env for AlphaGenome notebooks
│   └── r-motif-gsea-env.yml           # R env for complementary motif GSEA notebook
├── annotation/
│   └── gencode.v19.annotation.gtf.gz  # GENCODE v19 (hg19) — used by sashimi
├── notebooks/
│   ├── README.md
│   ├── 01_alphagenome_analysis/
│   │   ├── README.md
│   │   ├── alphagenome_cohort_analysis.ipynb          # cohort variant scoring
│   │   ├── alphagenome_genomAD_analysis.ipynb         # gnomAD panel scoring
│   │   ├── alphagenome_ISM_analysis.ipynb             # ISM on intron 6 hotspot (160 bp)
│   │   ├── COL4A5_ISM/                                # full-gene intron ISM workflow
│   │   │   ├── README.md                              # hg38 / RefSeq download + run order
│   │   │   ├── extract_introns.py                     # build intron BED/TSV/FASTA for NM_033380.3 from hg38
│   │   │   └── alphagenome_ISM_COL4A5.ipynb           # ISM on all COL4A5 introns
│   │   ├── alphagenome_intron6_cohort/                # cohort outputs (gitignored)
│   │   ├── alphagenome_intron6_genomAD/               # gnomAD outputs (gitignored)
│   │   ├── alphagenome_ISM_COL4A5_intron6/            # intron 6 ISM outputs (gitignored)
│   │   ├── alphagenome_ISM_COL4A5_intron47/           # intron 47 ISM outputs (gitignored)
│   │   └── alphagenome_ISM_COL4A5/                    # per-intron ISM outputs from COL4A5_ISM/ (gitignored)
│   └── 02_postprocessing_plots_summary/
│       ├── README.md
│       ├── alphagenome_analysis_notebook_ISM.ipynb            # intron 6 ISM summary + heatmaps
│       ├── alphagenome_analysis_notebook_ISM_intron47.ipynb   # intron 47 ISM summary + heatmaps
│       ├── alphagenome_analysis_notebook_cohort.ipynb         # cohort summary + heatmaps
│       ├── alphagenome_analysis_notebook_genomAD.ipynb        # gnomAD summary + heatmaps
│       ├── ISM_motif_annotation.ipynb                         # motif-level annotation
│       ├── motif_enrichment_pipeline.py                       # builds .rnk / .gmt / FASTA for GSEA
│       ├── run_complementary_motif_analysis.ipynb               # R-based fgsea/clusterProfiler cross-check
│       ├── plots_alphagenome_ISM/                               # intron 6 postprocessing outputs (gitignored)
│       ├── plots_alphagenome_ISM_intron47/                      # intron 47 postprocessing outputs (gitignored)
│       ├── plots_patients_cohort/                               # cohort postprocessing outputs (gitignored)
│       └── plots_alphagenome_genomAD/                           # gnomAD postprocessing outputs (gitignored)
└── targeted_rnaseq/                                             # ggsashimi inputs + figures
    ├── README.md
    ├── samples.tsv                                              # cohort 1 sample sheet (21 BAMs)
    ├── palette.tsv                                              # cohort 1 color palette
    ├── sashimi_COL4A5_intron6.{pdf,png}                           # cohort 1 sashimi figure
    ├── samples2.tsv                                               # cohort 2 sample sheet
    ├── palette2.tsv                                               # cohort 2 color palette
    ├── sashimi_COL4A5_intron6_inconclusive.pdf                  # cohort 2 sashimi figure
    ├── bams/                                                      # gitignored — see Zenodo (DOI below)
    └── bams2/                                                     # gitignored — cohort 2 BAMs (local only)
```

> **Targeted RNA-seq BAM files**: the 21 BAMs (+ `.bai` indexes) referenced by `targeted_rnaseq/samples.tsv` are archived on Zenodo: [10.5281/zenodo.19854065](https://doi.org/10.5281/zenodo.19854065). Download them into `targeted_rnaseq/bams/` (preserving filenames) before reproducing the sashimi plots.

## Workflow Overview

1. **`notebooks/01_alphagenome_analysis/`** — runs AlphaGenome on the cohort, gnomAD, and ISM variant sets and exports per-variant `csv/` and `png/` artifacts:
   - `alphagenome_cohort_analysis.ipynb` and `alphagenome_genomAD_analysis.ipynb` — patient cohort and gnomAD panels.
   - `alphagenome_ISM_analysis.ipynb` — saturation mutagenesis on the intron 6 hotspot (160 bp).
   - `COL4A5_ISM/` — download hg38 + UCSC RefSeq (`README.md`), extract all `NM_033380.3` introns with `extract_introns.py`, then run `alphagenome_ISM_COL4A5.ipynb` for per-intron ISM.
2. **`notebooks/02_postprocessing_plots_summary/`** — turns those CSVs into:
   - cohort/gnomAD/ISM summary tables, score heatmaps, and pathogenicity calls (`alphagenome_analysis_notebook_ISM.ipynb`, `alphagenome_analysis_notebook_ISM_intron47.ipynb`, `alphagenome_analysis_notebook_cohort.ipynb`, `alphagenome_analysis_notebook_genomAD.ipynb`);
   - motif-level inputs for GSEA (`motif_scores.rnk`, `motif_sets.gmt`, `motif_sequences.fasta`) via `motif_enrichment_pipeline.py`, with motif → variant traceability in `motif_scores_table.csv`;
   - an in-R GSEA cross-check via `run_complementary_motif_analysis.ipynb` (clusterProfiler/fgsea + enrichplot/ggseqlogo/pheatmap), output to `plots_alphagenome_ISM/complementary_motif_analysis/`.
3. **`targeted_rnaseq/`** — sample sheets, palettes, and sashimi figures for two RNA-seq cohorts (`samples.tsv` / `palette.tsv` and `samples2.tsv` / `palette2.tsv`). BAMs are gitignored due to size; cohort 1 BAMs are on Zenodo (see below).

## Conda Environments and JupyterLab
Two Conda environments are used:

```bash
# Python env for AlphaGenome notebooks and motif_enrichment_pipeline.py
conda env create -f coda/alphagenome-env.yml
conda activate alphagenome-env
python -m ipykernel install --user --name alphagenome-env --display-name "alphagenome-env"

# R env for the complementary motif GSEA notebook
conda env create -f coda/r-motif-gsea-env.yml
conda activate r-motif-gsea-env
# IRkernel registration is handled inside the notebook on first run

# Launch JupyterLab from the repo root
jupyter lab
```

Then open each notebook and select the matching kernel:
- Python notebooks → `alphagenome-env` (Python 3.11)
- `run_complementary_motif_analysis.ipynb` → `R (r-motif-gsea-env)`

## Reproducibility

Reproducibility details (inputs, versions, parameters, execution order, and outputs) are documented per stage:
- `CITATION.cff`
- `coda/README.md`
- `coda/alphagenome-env.yml`, `coda/r-motif-gsea-env.yml`
- `notebooks/README.md`
- `notebooks/01_alphagenome_analysis/README.md`
- `notebooks/01_alphagenome_analysis/COL4A5_ISM/README.md`
- `notebooks/02_postprocessing_plots_summary/README.md`
- `targeted_rnaseq/README.md` (sample design, hg19 reference, ggsashimi command)


