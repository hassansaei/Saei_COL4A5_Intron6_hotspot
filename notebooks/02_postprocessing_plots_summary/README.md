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

## Script in this folder

- `motif_enrichment_pipeline.py`
  - Input: `plots_alphagenome_ISM/ism_variant_summary.csv`
  - Outputs:
    - `motif_scores_table.csv`
    - `motif_scores.rnk`
    - `motif_sets.gmt`
    - `motif_sequences.fasta`
  - These files is used for GSEA preranked analysis and motif visualization.

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
7. (Optional) run GSEA preranked using the generated `.rnk` and `.gmt`, then review outputs in `motif_analysis.GseaPreranked.*`.

## Environment

- Use the Conda environment defined in `coda/alphagenome-env.yml`.
- Notebook kernel: `alphagenome-env` (Python 3.11).
