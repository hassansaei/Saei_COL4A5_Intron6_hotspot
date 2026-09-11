<p align="center">
  <img src="docs/logo.jpg" alt="Deep-intronic mutational hotspot schematic: exons as boxes, a concentrated intron hotspot, and a smaller secondary site" width="560">
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19854065.svg)](https://doi.org/10.5281/zenodo.19854065)

# Mapping deep intronic mutational hotspots by in silico mutagenesis enables single antisense oligonucleotide correction of multiple variants

**Authors:** Hassan Saei, Béatrice Ardin, Nicolas Kaiser, Mouad Ouahmane, Olivier Gribouval, Vincent Moriniere, Florian J. Wopperer, Korbinian M. Riedhammer, Daniel P. Gale, Omid Sadeghi-Alavijeh, Claire Goursaud, Olivier Grunewald, Chloe Prosper, Louis Lebreton, Marion Rabant, Carsten Bergmann, Corinne Antignac, Michael S. Wiesener, Guillaume Dorval

This repository contains scripts and notebooks used to generate figures for the manuscript on deep intronic variation hotspot(s) in `COL4A5` responsible for X-linked Alport syndrome. 

We organize the project for reproducibility, with separate stages for:
- AlphaGenome in-silico analyses (variant scoring, ISM saturation mutagenesis)
- Full-gene, gene-agnostic intron extraction and per-intron batch ISM (COL4A5, plus paralogs COL4A3/COL4A4)
- Downstream postprocessing, motif enrichment, and summary plotting
- Complementary in-R motif MSEA cross-check
- Targeted RNA-seq sashimi visualization of the intron 6 hotspot

## Glossary

| Term | Meaning in this repository |
|------|----------------------------|
| **AlphaGenome** | Deep-learning model for predicting variant effects on splicing and related molecular phenotypes; used here via the AlphaGenome API to score cohort, gnomAD, and ISM variants. See [References](#references) and [AlphaGenome docs](https://www.alphagenomedocs.com/). |
| **ISM** | In silico mutagenesis: the pipeline scores every single-nucleotide substitution in a region or intron in silico (saturation mutagenesis). |
| **MSEA** | Motif set enrichment analysis: preranked enrichment of 7-mer motif sets built from ISM splice-impact scores (`.rnk` + `.gmt` inputs). Implemented with GSEA-style preranked methods (`fgsea` / clusterProfiler). See [References](#references). |
| **ISS / ESE** | Intronic splicing silencer / exonic splicing enhancer: cis-regulatory motifs whose disruption can alter splicing. |
| **BH** | Benjamini–Hochberg multiple-testing correction that we apply to sliding-window hotspot tests within each intron. |
| **FDR** | False discovery rate; significance threshold for hotspot windows (default 0.05). |
| **NES** | Normalized enrichment score from MSEA output (direction and magnitude of motif-set enrichment). |
| **H/M** | High or Moderate impact ISM variants (splice-site max quantile ≥ 0.999 or ≥ 0.99). |

## References

Methods and software cited in this repository:

- **AlphaGenome** (variant scoring and ISM): Avsec, Ž., Latysheva, N., Cheng, J. et al. Advancing regulatory variant effect prediction with AlphaGenome. *Nature* **649**, 1206–1218 (2026). https://doi.org/10.1038/s41586-025-10014-0  
  API and Python SDK: https://www.alphagenomedocs.com/

- **GSEA** (gene set enrichment analysis; preranked enrichment used for MSEA): Subramanian, A., Tamayo, P., Mootha, V.K. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc. Natl Acad Sci. USA* **102**, 15545–15550 (2005). https://doi.org/10.1073/pnas.0506580102  
  In this repo, preranked MSEA runs use `fgsea` (Korotkevich, G. et al., *bioRxiv* 2019, https://doi.org/10.1101/060012) and `clusterProfiler` (Wu, T. et al., *Innovation* **2**, 100141, 2021, https://doi.org/10.1016/j.xinn.2021.100141) via `run_complementary_motif_analysis.ipynb`.

## Data availability
We archived targeted RNA-seq BAM files (+ `.bai` indexes) on Zenodo:

**DOI:** [10.5281/zenodo.19854065](https://doi.org/10.5281/zenodo.19854065)

Download the archive and place the files in `targeted_rnaseq/bams/` (preserving filenames) before reproducing the sashimi figure. See `targeted_rnaseq/README.md` for sample design and ggsashimi commands.

## Repository Structure

```text
.
├── README.md
├── CITATION.cff
├── LICENSE
├── docs/
│   └── logo.jpg                       # repository logo (intron hotspot schematic)
├── coda/                              # Conda environment specs
│   ├── README.md
│   ├── alphagenome-env.yml            # Python env for AlphaGenome notebooks
│   └── r-motif-gsea-env.yml           # R env for complementary motif MSEA notebook
├── annotation/
│   └── gencode.v19.annotation.gtf.gz  # GENCODE v19 (hg19) — used by sashimi
├── notebooks/
│   ├── README.md
│   ├── 01_alphagenome_analysis/
│   │   ├── README.md
│   │   ├── alphagenome_cohort_analysis.ipynb          # cohort variant scoring
│   │   ├── alphagenome_genomAD_analysis.ipynb         # gnomAD panel scoring
│   │   ├── alphagenome_ISM_analysis.ipynb             # ISM on intron 6 hotspot (160 nt)
│   │   ├── ISM/                                       # full-gene, gene-agnostic intron ISM workflow
│   │   │   ├── README.md                              # RefSeq/hg38 download + CLI + run order
│   │   │   ├── extract_introns.py                     # CLI: build intron BED/TSV/FASTA for any RefSeq transcript
│   │   │   └── alphagenome_ISM_batch.ipynb           # batch VCF ISM on all introns of a gene (COL4A5/COL4A3/COL4A4)
│   │   ├── alphagenome_intron6_cohort/                # cohort outputs (gitignored)
│   │   ├── alphagenome_intron6_genomAD/               # gnomAD outputs (gitignored)
│   │   ├── alphagenome_ISM_COL4A5_intron6/            # intron 6 ISM outputs (gitignored)
│   │   ├── alphagenome_ISM_COL4A5_intron47/           # intron 47 ISM outputs (gitignored)
│   │   ├── alphagenome_ISM_COL4A5/                    # per-intron ISM outputs for COL4A5 (gitignored)
│   │   ├── alphagenome_ISM_COL4A3/                    # per-intron ISM outputs for COL4A3 (gitignored)
│   │   └── alphagenome_ISM_COL4A4/                    # per-intron ISM outputs for COL4A4 (gitignored)
│   └── 02_postprocessing_plots_summary/
│       ├── README.md
│       ├── alphagenome_analysis_notebook_ISM.ipynb            # intron 6 ISM summary + heatmaps
│       ├── alphagenome_analysis_notebook_ISM_intron47.ipynb   # intron 47 ISM summary + heatmaps
│       ├── alphagenome_analysis_notebook_ISM_all_introns.ipynb # all-introns ISM summary + cross-intron hotspot comparison
│       ├── alphagenome_analysis_notebook_cohort.ipynb         # cohort summary + heatmaps
│       ├── alphagenome_analysis_notebook_genomAD.ipynb        # gnomAD summary + heatmaps
│       ├── ISM_motif_annotation.ipynb                         # motif-level annotation
│       ├── motif_enrichment_pipeline.py                       # builds .rnk / .gmt / FASTA for MSEA (per hotspot intron)
│       ├── plot_gsea_cohort_summary.py                        # cohort MSEA NES heatmap across hotspot introns
│       ├── plot_hotspot_motif_disruptions.py                  # strict ISS/ESE motif-disruption counts per hotspot
│       ├── run_complementary_motif_analysis.ipynb               # R-based fgsea/clusterProfiler cross-check
│       ├── plots_alphagenome_ISM/                               # legacy intron 6 single-window postprocessing (gitignored)
│       ├── plots_alphagenome_ISM_Intron<N>/                     # per-hotspot-intron motif/MSEA outputs (Intron6, etc.; gitignored)
│       ├── plots_alphagenome_ISM_intron47/                      # intron 47 postprocessing outputs (gitignored)
│       ├── plots_alphagenome_ISM_<GENE>_all_introns/            # all-introns postprocessing outputs, per gene (gitignored)
│       ├── motif_enrichment_cohort_summary/                     # cohort MSEA summary plots (gitignored)
│       ├── motif_disruption_cohort_summary/                     # cohort motif-disruption plots (gitignored)
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

> **Targeted RNA-seq BAM files**: the 21 BAMs (+ `.bai` indexes) referenced by `targeted_rnaseq/samples.tsv` sit on Zenodo: [10.5281/zenodo.19854065](https://doi.org/10.5281/zenodo.19854065). Download them into `targeted_rnaseq/bams/` (preserving filenames) before reproducing the sashimi plots.

## Workflow Overview

1. **`notebooks/01_alphagenome_analysis/`** — runs AlphaGenome on the cohort, gnomAD, and ISM variant sets and exports per-variant `csv/` and `png/` artifacts:
   - `alphagenome_cohort_analysis.ipynb` and `alphagenome_genomAD_analysis.ipynb` — patient cohort and gnomAD panels.
   - `alphagenome_ISM_analysis.ipynb` — saturation mutagenesis on the intron 6 hotspot (160 nt).
   - `ISM/` — gene-agnostic full-gene workflow: download hg38 + UCSC RefSeq (`README.md`), extract all introns of a transcript with the `extract_introns.py` CLI (e.g. `python extract_introns.py NM_033380.3`), then run `alphagenome_ISM_batch.ipynb` for fast batch VCF-based per-intron ISM (defaults to COL4A5; also run for COL4A3/COL4A4).
2. **`notebooks/02_postprocessing_plots_summary/`** — turns those CSVs into:
   - cohort/gnomAD/ISM summary tables, score heatmaps, and pathogenicity calls (`alphagenome_analysis_notebook_ISM.ipynb`, `alphagenome_analysis_notebook_ISM_intron47.ipynb`, `alphagenome_analysis_notebook_cohort.ipynb`, `alphagenome_analysis_notebook_genomAD.ipynb`);
   - a gene-agnostic all-introns summary with cross-intron hotspot comparison (`alphagenome_analysis_notebook_ISM_all_introns.ipynb`, consuming a `../01_alphagenome_analysis/alphagenome_ISM_<GENE>/` batch run); hotspot **concentration** uses one representative interior hotspot per intron (`CONCENTRATION_MODE = 'top_hotspot'`);
   - per-hotspot-intron motif inputs for MSEA (`motif_scores.rnk`, `motif_sets.gmt`, `motif_sequences.fasta`) via `motif_enrichment_pipeline.py` (hotspot introns 1, 4, 6, 30, 44, 49), with motif-to-variant traceability in `motif_scores_table.csv`;
   - an in-R MSEA cross-check per hotspot intron via `run_complementary_motif_analysis.ipynb` (clusterProfiler/fgsea + enrichplot/ggseqlogo/pheatmap), output to `plots_alphagenome_ISM_Intron<N>/complementary_motif_analysis/`;
   - cohort-level summary figures via `plot_gsea_cohort_summary.py` (MSEA NES heatmap) and `plot_hotspot_motif_disruptions.py` (strict ISS/ESE motif-disruption counts).
3. **`targeted_rnaseq/`** — sample sheets, palettes, and sashimi figures for two RNA-seq cohorts (`samples.tsv` / `palette.tsv` and `samples2.tsv` / `palette2.tsv`). We gitignore BAMs due to size; cohort 1 BAMs are on Zenodo (see below).

## Conda Environments and JupyterLab
We use two Conda environments:

```bash
# Python env for AlphaGenome notebooks and motif_enrichment_pipeline.py
conda env create -f coda/alphagenome-env.yml
conda activate alphagenome-env
python -m ipykernel install --user --name alphagenome-env --display-name "alphagenome-env"

# R env for the complementary motif MSEA notebook
conda env create -f coda/r-motif-gsea-env.yml
conda activate r-motif-gsea-env
# The notebook registers the IRkernel on first run

# Launch JupyterLab from the repo root
jupyter lab
```

Then open each notebook and select the matching kernel:
- Python notebooks: use `alphagenome-env` (Python 3.11)
- `run_complementary_motif_analysis.ipynb`: use `R (r-motif-gsea-env)`

## Reproducibility

We document reproducibility details (inputs, versions, parameters, execution order, and outputs) per stage in:
- `CITATION.cff`
- `coda/README.md`
- `coda/alphagenome-env.yml`, `coda/r-motif-gsea-env.yml`
- `notebooks/README.md`
- `notebooks/01_alphagenome_analysis/README.md`
- `notebooks/01_alphagenome_analysis/ISM/README.md`
- `notebooks/02_postprocessing_plots_summary/README.md`
- `targeted_rnaseq/README.md` (sample design, hg19 reference, ggsashimi command)


