# Notebooks

This directory contains all analysis notebooks used for manuscript figures, split by workflow stage.

## Folder Organization

- `01_alphagenome_analysis/`  
  Notebooks that run AlphaGenome analyses and produce primary result files
  (cohort, gnomAD, intron 6 hotspot ISM, and the gene-agnostic full-gene batch
  ISM workflow in `ISM/`). See `01_alphagenome_analysis/README.md` and
  `01_alphagenome_analysis/ISM/README.md`.

- `02_postprocessing_plots_summary/`  
  Notebooks that load AlphaGenome outputs and generate processed tables, plots,
  summary statistics, and motif-enrichment inputs for figures. See
  `02_postprocessing_plots_summary/README.md`.

## Suggested Execution Flow

1. Run notebooks in `01_alphagenome_analysis/`
2. Verify generated outputs
3. Run notebooks in `02_postprocessing_plots_summary/`

## Reproducibility Checklist

- Project title: COL4A5 intron 6 pseudoexon hotspot
- Manuscript title: Recurrent pseudoexon activation in COL4A5 intron 6 defines a therapeutically actionable hotspot
- Corresponding author/contact: Hassan Saei / Guillaume Dorval (see root `README.md` and `CITATION.cff`)
- Date last updated: 2026-07-05

### Environment

- Python version: 3.11 (kernel `alphagenome-env`); R notebook uses kernel `R (r-motif-gsea-env)`
- Package manager: Conda
- Environment file path: `coda/alphagenome-env.yml` (Python) and `coda/r-motif-gsea-env.yml` (R)
- OS tested: macOS

### Inputs

- Reference genome build: hg38 (GRCh38); GENCODE v46 annotations for AlphaGenome
- Gene/transcript definitions: COL4A5 `NM_033380.3` (plus paralogs COL4A3 `NM_000091.5`, COL4A4 `NM_000092.5`) via the UCSC `ncbiRefSeq` table; introns extracted with `01_alphagenome_analysis/ISM/extract_introns.py`
- Variant list source: patient cohort, gnomAD panel, and in-silico saturation mutagenesis (all single-nucleotide substitutions per region/intron)
- Any preprocessing applied: tissue filtering to kidney (`UBERON:0002113`); per-intron VCFs built for batch scoring

### Outputs

- Primary AlphaGenome output location: `01_alphagenome_analysis/alphagenome_intron6_cohort/`, `alphagenome_intron6_genomAD/`, `alphagenome_ISM_COL4A5_intron6/`, `alphagenome_ISM_COL4A5_intron47/`, and per-gene `alphagenome_ISM_<GENE>/` (all gitignored)
- Processed table output location: `02_postprocessing_plots_summary/plots_*/` (e.g. `plots_alphagenome_ISM/`, `plots_alphagenome_ISM_<GENE>_all_introns/`, gitignored)
- Figure output location: the same `plots_*/` folders (PNG/SVG)

### Notebook Run Order

1. `01_alphagenome_analysis/alphagenome_cohort_analysis.ipynb`
2. `01_alphagenome_analysis/alphagenome_genomAD_analysis.ipynb`
3. `01_alphagenome_analysis/alphagenome_ISM_analysis.ipynb` (intron 6 hotspot)
4. `01_alphagenome_analysis/ISM/alphagenome_ISM_batch.ipynb` (full-gene batch ISM; run after `extract_introns.py`)
5. `02_postprocessing_plots_summary/alphagenome_analysis_notebook_cohort.ipynb`
6. `02_postprocessing_plots_summary/alphagenome_analysis_notebook_genomAD.ipynb`
7. `02_postprocessing_plots_summary/alphagenome_analysis_notebook_ISM.ipynb` (intron 6) and `alphagenome_analysis_notebook_ISM_intron47.ipynb` (intron 47)
8. `02_postprocessing_plots_summary/alphagenome_analysis_notebook_ISM_all_introns.ipynb` (all introns + cross-intron hotspot comparison)
9. `02_postprocessing_plots_summary/motif_enrichment_pipeline.py`, then `ISM_motif_annotation.ipynb` and (optional) `run_complementary_motif_analysis.ipynb`
