# Recurrent pseudoexon activation in COL4A5 intron 6 defines a therapeutically actionable hotspot

This repository contains scripts and notebooks used to generate figures for the manuscript on deep intronic variation hotspot(s) in `COL4A5` responsible for X-linked Alport syndrome. 

The project is organized for reproducibility and clear separation between:
- AlphaGenome in-silico analyses (variant scoring, ISM saturation mutagenesis)
- Downstream postprocessing, motif enrichment, and summary plotting
- Complementary in-R motif GSEA cross-check
- Targeted RNA-seq sashimi visualization of the intron 6 hotspot

## Repository Structure

```text
.
├── README.md
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
│   │   ├── alphagenome_cohort_analysis.ipynb     # cohort variant scoring
│   │   ├── alphagenome_genomAD_analysis.ipynb    # gnomAD panel scoring
│   │   └── alphagenome_ISM_analysis.ipynb        # in-silico saturation mutagenesis
│   └── 02_postprocessing_plots_summary/
│       ├── README.md
│       ├── alphagenome_analysis_notebook_ISM.ipynb       # ISM summary + heatmaps
│       ├── alphagenome_analysis_notebook_cohort.ipynb    # cohort summary + heatmaps
│       ├── alphagenome_analysis_notebook_genomAD.ipynb   # gnomAD summary + heatmaps
│       ├── ISM_motif_annotation.ipynb                    # motif-level annotation
│       ├── motif_enrichment_pipeline.py                  # builds .rnk / .gmt / FASTA for GSEA
│       └── run_complementary_motif_analysis.ipynb        # R-based fgsea/clusterProfiler cross-check
└── Targeted_rnaseq/                                      # ggsashimi inputs + figure
    ├── README.md
    ├── samples.tsv                                       # ggsashimi -b sample sheet (21 BAMs)
    ├── palette.tsv                                       # ggsashimi -P color palette
    ├── sashimi_COL4A5_intron6.{pdf,png}                  # rendered sashimi figure
    └── bams/                                             # gitignored — see Zenodo (DOI below)
```

> **Targeted RNA-seq BAM files**: the 21 BAMs (+ `.bai` indexes) referenced by `Targeted_rnaseq/samples.tsv` are archived on Zenodo: [10.5281/zenodo.19854065](https://doi.org/10.5281/zenodo.19854065). Download them into `Targeted_rnaseq/bams/` (preserving filenames) before reproducing the sashimi figure.

## Workflow Overview

1. **`notebooks/01_alphagenome_analysis/`** — runs AlphaGenome on the cohort, gnomAD, and exhaustive ISM variant sets and exports per-variant `csv/` and `png/` artifacts.
2. **`notebooks/02_postprocessing_plots_summary/`** — turns those CSVs into:
   - cohort/gnomAD/ISM summary tables, score heatmaps, and pathogenicity calls;
   - motif-level inputs for GSEA (`motif_scores.rnk`, `motif_sets.gmt`, `motif_sequences.fasta`) via `motif_enrichment_pipeline.py`, with motif → variant traceability in `motif_scores_table.csv`;
   - an in-R GSEA cross-check via `run_complementary_motif_analysis.ipynb` (clusterProfiler/fgsea + enrichplot/ggseqlogo/pheatmap), output to `plots_alphagenome_ISM/complementary_motif_analysis/`.
3. **`Targeted_rnaseq/`** — sample sheet, palette, and the `sashimi_COL4A5_intron6` figure produced from BAMs (BAMs themselves are gitignored due to size).

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
- `coda/README.md`
- `coda/alphagenome-env.yml`, `coda/r-motif-gsea-env.yml`
- `notebooks/README.md`
- `notebooks/01_alphagenome_analysis/README.md`
- `notebooks/02_postprocessing_plots_summary/README.md`
- `Targeted_rnaseq/README.md` (sample design, hg19 reference, ggsashimi command)


