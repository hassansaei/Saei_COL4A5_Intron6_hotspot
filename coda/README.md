# Conda Environment Setup (`alphagenome-env`)

This folder contains a reproducible Conda environment file for the AlphaGenome notebook workflow.

## Files

- `alphagenome-env.yml`: Conda environment definition used by the notebook kernel `alphagenome-env`.

## 1) Create the environment

From the repository root:

```bash
conda env create -f coda/alphagenome-env.yml
```

## 2) Activate the environment

```bash
conda activate alphagenome-env
```

## 3) Register the Jupyter kernel

```bash
python -m ipykernel install --user --name alphagenome-env --display-name "alphagenome-env"
```

## 4) Launch Jupyter

```bash
jupyter lab
```

Open notebooks and select kernel: `alphagenome-env`.

## 5) Verify key packages

```bash
python -c "import alphagenome, pandas, numpy, matplotlib, tqdm; print('Environment OK')"
```

## Update the environment (after changing YAML)

```bash
conda env update -n alphagenome-env -f coda/alphagenome-env.yml --prune
```

## Remove the environment (optional)

```bash
conda remove -n alphagenome-env --all
```
