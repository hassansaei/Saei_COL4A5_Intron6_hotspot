# AlphaGenome Analysis Notebooks

This folder contains the primary AlphaGenome workflows for the `COL4A5` intron 6 hotspot analyses.

## Notebooks in this folder

- `alphagenome_cohort_analysis.ipynb`
  - Cohort-focused variant analysis (19 variants).
  - Exports per-variant plots and per-variant splice scorer CSV files.
- `alphagenome_genomAD_analysis.ipynb`
  - gnomAD variant panel analysis (14 variants).
  - Exports per-variant plots and per-variant splice scorer CSV files.
- `alphagenome_ISM_analysis.ipynb`
  - Comprehensive saturation mutagenesis (all substitutions across input sequence).
  - Exports per-variant plots and per-variant splice scorer CSV files.

See `ISM/README.md` for the full-gene, gene-agnostic intron ISM workflow
(`extract_introns.py` CLI + `alphagenome_ISM_batch.ipynb`), which builds a
per-intron VCF of all single-nucleotide substitutions and batch-scores them.
Configure `GENE_NAME` / `TRANSCRIPT_ID` / `ONTOLOGY_CURIE` in the notebook
(defaults COL4A5 / NM_033380.3 / kidney; also run for COL4A3 and COL4A4). Outputs
land in per-gene folders `alphagenome_ISM_<GENE>/`.

## Current output/result folders

- `alphagenome_intron6_cohort/`
  - `png/`
  - `csv/`
  - `variant_run_summary.csv`
- `alphagenome_intron6_genomAD/`
  - `png/`
  - `csv/`
  - `variant_run_summary.csv`
- `alphagenome_ISM_COL4A5_intron6/`
  - `png/`
  - `csv/`
  - (summary table is returned in notebook; save to CSV if needed)
- `alphagenome_ISM_<GENE>/` (e.g. `alphagenome_ISM_COL4A5/`, `alphagenome_ISM_COL4A3/`, `alphagenome_ISM_COL4A4/`)
  - one `intron_*/` subfolder per intron, each with `<intron>_variants.vcf` and `csv/alphagenome_scores_<intron>.csv`
  - produced by `ISM/alphagenome_ISM_batch.ipynb`

## Output naming convention

- Full-window plot (4500 bp): `alphagenome_variant_<variant>_kidney.png`
- Zoomed plot (1000 bp): `alphagenome_variant_<variant>_kidney_zoomed.png`
- Variant scorer table: `alphagenome_scores_<variant>.csv`
- Run summary table (cohort/gnomAD): `variant_run_summary.csv`

Variant format follows:
- `g.<position><REF>><ALT>`
- examples: `g.108570705G>C`, `g.108570778G>GC`

## Reproducibility notes

- Use the Conda environment in `coda/alphagenome-env.yml`.
- Kernel used by notebooks: `alphagenome-env` (Python 3.11).
- Reference annotation source in notebooks:
  - GENCODE v46 hg38 feather from Google Cloud Storage.
- Ontology term used for tissue filtering in these analyses:
  - `UBERON:0002113` (kidney).

## Recommended run order

1. Run setup cells (imports, API client, GTF loading, transcript extractors).
2. Run `alphagenome_cohort_analysis.ipynb`.
3. Run `alphagenome_genomAD_analysis.ipynb`.
4. Run `alphagenome_ISM_analysis.ipynb`.
5. Save/verify summary outputs, then continue to notebooks in `../02_postprocessing_plots_summary/`.

## Practical notes

- In long batch runs, Matplotlib may warn about many open figures.
- If memory usage grows, add `plt.close(plot)` and `plt.close(plot_zoomed)` after each `savefig` in batch loops.
