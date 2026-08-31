# Conda Environments

This folder contains the reproducible Conda environment files used across
the project. Two separate environments are maintained so the Python
AlphaGenome stack and the R motif-MSEA stack stay isolated.

Term definitions (including **MSEA**, motif set enrichment analysis): see the
[Glossary](../README.md#glossary) in the repository root README.

## Files

- `alphagenome-env.yml` — Python env (Python 3.11) for the AlphaGenome
  notebooks in `notebooks/01_alphagenome_analysis/` and the postprocessing
  notebooks/scripts in `notebooks/02_postprocessing_plots_summary/`
  (including `motif_enrichment_pipeline.py`). Notebook kernel:
  `alphagenome-env`.
- `r-motif-gsea-env.yml` — R env (R ≥ 4.3) for the complementary motif
  MSEA notebook
  `notebooks/02_postprocessing_plots_summary/run_complementary_motif_analysis.ipynb`.
  Provides `r-base`, `r-irkernel`, `r-biocmanager`, `r-ggplot2`,
  `r-pheatmap`, `r-ggseqlogo`, `r-gdtools`, plus build toolchain
  (`c-/cxx-/fortran-compiler`, `make`, `pkg-config`) needed to compile
  Bioconductor packages installed at first run (`clusterProfiler`,
  `enrichplot`, `fgsea`). Notebook kernel: `R (r-motif-gsea-env)`.

> **Where each notebook expects its kernel** is documented in
> `notebooks/02_postprocessing_plots_summary/README.md`.

---

## A. Python env — `alphagenome-env`

### 1) Create

From the repository root:

```bash
conda env create -f coda/alphagenome-env.yml
```

### 2) Activate

```bash
conda activate alphagenome-env
```

### 3) Register the Jupyter kernel

```bash
python -m ipykernel install --user --name alphagenome-env --display-name "alphagenome-env"
```

### 4) Launch Jupyter

```bash
jupyter lab
```

Open notebooks and select kernel: `alphagenome-env`.

### 5) Verify key packages

```bash
python -c "import alphagenome, pandas, numpy, matplotlib, tqdm; print('alphagenome-env OK')"
```

`motif_enrichment_pipeline.py` additionally requires `biopython`
(`from Bio.Seq import Seq`); install via
`pip install biopython` inside the env if it is not already present.

### 6) Update / remove

```bash
conda env update -n alphagenome-env -f coda/alphagenome-env.yml --prune
conda remove   -n alphagenome-env --all
```

---

## B. R env — `r-motif-gsea-env`

### 1) Create and activate

```bash
conda env create -f coda/r-motif-gsea-env.yml
conda activate r-motif-gsea-env
```

### 2) Register the IRkernel

The first cells of `run_complementary_motif_analysis.ipynb` install any
missing Bioconductor packages (`clusterProfiler`, `enrichplot`, `fgsea`)
and register the IRkernel automatically. To register manually instead:

```bash
R -e 'IRkernel::installspec(name = "r-motif-gsea-env", displayname = "R (r-motif-gsea-env)")'
```

### 3) Launch Jupyter and select the kernel

```bash
jupyter lab
```

Open `run_complementary_motif_analysis.ipynb` and select kernel
`R (r-motif-gsea-env)`.

### 4) Verify key packages

```bash
R -e 'suppressPackageStartupMessages({library(clusterProfiler); library(fgsea); library(enrichplot); library(ggseqlogo); library(pheatmap)}); cat("r-motif-gsea-env OK\n")'
```

### 5) Update / remove

```bash
conda env update -n r-motif-gsea-env -f coda/r-motif-gsea-env.yml --prune
conda remove   -n r-motif-gsea-env --all
```
