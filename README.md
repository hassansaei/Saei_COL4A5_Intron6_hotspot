# A deep intronic mutational hotspot in COL4A5 enables RNA-based personalized therapy in Alport syndrome

This repository contains scripts and notebooks used to generate figures for the manuscript on deep intronic variation hotspot(s) in `COL4A5` responsible for X-linked Alport syndrome. 

The project is organized for reproducibility and clear separation between:
- AlphaGenome analysis notebooks
- Downstream result processing, plotting, and summary notebooks

## Repository Structure

```text
.
├── coda/
│   ├── README.md
│   └── alphagenome-env.yml
├── LICENSE
├── README.md
└── notebooks/
    ├── README.md
    ├── 01_alphagenome_analysis/
    │   └── README.md
    └── 02_postprocessing_plots_summary/
        └── README.md
```

## Conda Environment and JupyterLab

Create and use the `alphagenome-env` environment from the repository root.

```bash
# Create environment and activate it
conda env create -f coda/alphagenome-env.yml
conda activate alphagenome-env
# Register kernel
python -m ipykernel install --user --name alphagenome-env --display-name "alphagenome-env"
# Launch JupyterLab
jupyter lab
```

Then open your notebook and select kernel `alphagenome-env`.

## Reproducibility

Reproducibility details (inputs, versions, parameters, execution order, and outputs) should be documented in:
- `coda/README.md`
- `coda/alphagenome-env.yml`
- `notebooks/README.md`
- `notebooks/01_alphagenome_analysis/README.md`
- `notebooks/02_postprocessing_plots_summary/README.md`


