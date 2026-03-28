# Notebooks

This directory contains all analysis notebooks used for manuscript figures, split by workflow stage.

## Folder Organization

- `01_alphagenome_analysis/`  
  Notebooks that run AlphaGenome analyses and produce primary result files.

- `02_postprocessing_plots_summary/`  
  Notebooks that load AlphaGenome outputs and generate processed tables, plots, and summary statistics for figures.

## Suggested Execution Flow

1. Run notebooks in `01_alphagenome_analysis/`
2. Verify generated outputs
3. Run notebooks in `02_postprocessing_plots_summary/`

## Reproducibility Checklist (fill in)

- Project title:
- Manuscript title:
- Corresponding author/contact:
- Date last updated:

### Environment

- Python version:
- Package manager:
- Environment file path (`environment.yml` or `requirements.txt`):
- OS tested:

### Inputs

- Reference genome build:
- Gene/transcript definitions:
- Variant list source:
- Any preprocessing applied:

### Outputs

- Primary AlphaGenome output location:
- Processed table output location:
- Figure output location:

### Notebook Run Order

- `01_alphagenome_analysis/<notebook_name>.ipynb`
- `02_postprocessing_plots_summary/<notebook_name>.ipynb`

Add additional notebook paths as needed.
