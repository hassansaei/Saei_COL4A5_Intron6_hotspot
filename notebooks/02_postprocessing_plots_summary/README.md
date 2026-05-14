# Postprocessing, Plots, and Summary Notebooks

This folder contains downstream postprocessing notebooks for the `COL4A5` intron 6 hotspot project.  
It turns primary AlphaGenome CSV outputs from `../01_alphagenome_analysis/` into summary tables, heatmaps, and motif-analysis inputs.

## Notebooks in this folder

- `alphagenome_analysis_notebook_ISM.ipynb`
  - Processes ISM per-variant score CSVs.
  - Builds variant-level summary table and heatmaps.
  - Current heatmap logic is configured to plot top variants in `plot_scores_heatmap`.
- `alphagenome_analysis_notebook_cohort.ipynb`
  - Processes patient cohort AlphaGenome score outputs.
  - Produces score heatmap, genomic-position plot, and pathogenicity summary CSV.
- `alphagenome_analysis_notebook_genomAD.ipynb`
  - Processes gnomAD variant panel AlphaGenome score outputs.
  - Produces score heatmap, genomic-position plot, and pathogenicity summary CSV.
- `ISM_motif_annotation.ipynb`
  - Annotates motif-level ISM outputs.
  - Reads `motif_scores_table.csv` and writes `motif_scores_annotated.csv`.
- `run_complementary_motif_analysis.ipynb` (R notebook, kernel `R (r-motif-gsea-env)`)
  - Independent in-R reimplementation of motif GSEA, complementary to the Python/CLI GSEA preranked run.
  - Reads `plots_alphagenome_ISM/motif_scores.rnk` and `plots_alphagenome_ISM/motif_sets.gmt`.
  - Uses `clusterProfiler::GSEA` (via `fgsea`), with `enrichplot`, `ggseqlogo`, and `pheatmap` for visualization.
  - Writes all outputs to `plots_alphagenome_ISM/complementary_motif_analysis/`:
    - `gsea_results_r.csv`, `gsea_significant_sets_r.csv`
    - `gsea_top_enriched_sets_r.csv`, `gsea_top_depleted_sets_r.csv`
    - `gsea_top_sets_nes_r.png`
    - `gsea_curve_enriched_<set>_r.png` and `gsea_curve_depleted_<set>_r.png` (top-N curves, default N = 5)
    - `leading_edge_membership_heatmap_r.png` and a leading-edge sequence-logo PNG
  - Tunable parameters in the first cell: `min_size`, `max_size`, `top_curves`, RNG seed.

## Script in this folder

- `motif_enrichment_pipeline.py`
  - Input: `plots_alphagenome_ISM/ism_variant_summary.csv`
  - Builds a 7-mer (`Context7mer`) variant context using the hard-coded intron-6 `REFERENCE_SEQ` and `region_start = 108570633`; window size is controlled by `WINDOW` (5–8).
  - Ranks each 7-mer by the **mean of per-variant `Splice Site Max Quantile`** among variants that create that context. The exported column is **`Mean Splice Site Max Quantile`** (not “mean mean”: it is one mean, over variants, of splice-site max quantiles). Configure source/label columns in `motif_enrichment_pipeline.py` via `VARIANT_SCORE_COLUMN` and `AGG_SCORE_COLUMN`.
  - Outputs (written next to the script, in the working directory):
    - `motif_scores_table.csv` — now includes a `Creating Variants` column listing the unique `Variant ID`s that produced each motif (preserves motif → variant traceability).
    - `motif_scores.rnk` — two-column GSEA preranked input (`motif`, score).
    - `motif_sets.gmt` — motif sets built from non-positional 3-mers, simple homopolymer repeats (`X-repeat`), and `CG-core`/`AG-core` families.
    - `motif_sequences.fasta` — one record per unique motif, used by MEME / `ggseqlogo`.
  - These files feed both the CLI GSEA preranked workflow and the R notebook above.

## Input dependencies

Expected upstream inputs come from `../01_alphagenome_analysis/`:

- ISM input folder:
  - `../01_alphagenome_analysis/alphagenome_ISM_COL4A5_intron6/csv/`
- Cohort input folder:
  - `../01_alphagenome_analysis/alphagenome_intron6_cohort/csv/`
- gnomAD input folder:
  - `../01_alphagenome_analysis/alphagenome_intron6_genomAD/csv/`

## Output folders and current artifacts

- `plots_alphagenome_ISM/`
  - `ism_variant_summary.csv`
  - `ism_heatmap.svg` (and notebook also saves PNG when run)
  - `variant_scores_heatmap.svg` (and notebook also saves PNG when run)
  - `motif_scores_table.csv`
  - `motif_scores_annotated.csv`
  - `motif_scores.rnk`
  - `motif_sets.gmt`
  - `motif_sequences.fasta`
  - `motif_sequence_high_impact.fasta`
  - `gsea_motif_analysis/` — per-motif-set HTML/TSV reports from the CLI GSEA preranked run (one `<SET>.html` + `<SET>.tsv` per gene set, e.g. `CG-CORE.*`, `ACG.*`).
  - `complementary_motif_analysis/` — outputs of `run_complementary_motif_analysis.ipynb` (R-based GSEA + leading-edge heatmap + seqlogo); see notebook section above for the full file list.
- `plots_patients_cohort/`
  - `variant_scores_heatmap.svg`
  - `variant_genomic_positions.svg`
  - `variant_pathogenicity_summary.csv`
- `plots_alphagenome_genomAD/`
  - `variant_scores_heatmap.svg`
  - `variant_genomic_positions.svg`
  - `variant_pathogenicity_summary.csv`

## Motif enrichment output

- `motif_analysis.GseaPreranked.<timestamp>/`
  - Contains GSEA report pages, `edb/` files, and enrichment tables generated from preranked motif analysis.

## Recommended run order

1. Run upstream analyses in `../01_alphagenome_analysis/`.
2. Run `alphagenome_analysis_notebook_ISM.ipynb`.
3. Run `alphagenome_analysis_notebook_cohort.ipynb`.
4. Run `alphagenome_analysis_notebook_genomAD.ipynb`.
5. Run `motif_enrichment_pipeline.py`.
6. Run `ISM_motif_annotation.ipynb`.
7. (Optional) run GSEA preranked using the generated `.rnk` and `.gmt`, then review outputs in `motif_analysis.GseaPreranked.*` and `plots_alphagenome_ISM/gsea_motif_analysis/`.
8. (Optional) Run `run_complementary_motif_analysis.ipynb` for an in-R GSEA cross-check; outputs land in `plots_alphagenome_ISM/complementary_motif_analysis/`.

## Environment

- Python notebooks/script: use the Conda environment defined in `coda/alphagenome-env.yml`. Notebook kernel: `alphagenome-env` (Python 3.11).
- `motif_enrichment_pipeline.py` additionally requires `biopython` (`from Bio.Seq import Seq`).
- `run_complementary_motif_analysis.ipynb` runs on a separate R kernel `R (r-motif-gsea-env)` with Bioconductor packages `clusterProfiler`, `enrichplot`, `fgsea`, and CRAN packages `ggplot2`, `ggseqlogo`, `pheatmap`, `gridExtra` (install commands are commented in the first setup cell).
