# Postprocessing, Plots, and Summary Notebooks

This folder contains downstream postprocessing notebooks and scripts for the `COL4A5` deep-intronic hotspot project.  
It turns primary AlphaGenome CSV outputs from `../01_alphagenome_analysis/` into summary tables, heatmaps, motif-analysis inputs, and cohort-level figures.

The workflow spans:
- **Single-window ISM** (intron 6 hotspot, 160 nt; intron 47 control)
- **Full-gene batch ISM** across all COL4A5 introns, with cross-intron hotspot comparison
- **Per-hotspot-intron MSEA** (motif set enrichment analysis; introns 1, 4, 6, 30, 44, 49) and cohort summaries

Term definitions (AlphaGenome, ISM, MSEA, ISS/ESE, BH, FDR, NES, H/M): see the [Glossary](../../README.md#glossary) in the repository root README.

## Notebooks in this folder

- `alphagenome_analysis_notebook_ISM.ipynb`
  - Processes ISM per-variant score CSVs (COL4A5 **intron 6**, 160 nt).
  - Builds variant-level summary table and heatmaps.
  - Current heatmap logic is configured to plot top variants in `plot_scores_heatmap`.
  - Output: `plots_alphagenome_ISM/` (legacy single-window path; motif/MSEA for intron 6 uses `plots_alphagenome_ISM_Intron6/`)
- `alphagenome_analysis_notebook_ISM_intron47.ipynb`
  - Same workflow for COL4A5 **intron 47** (123 nt; `chrX:108683908–108684030`).
  - Input: `../01_alphagenome_analysis/alphagenome_ISM_COL4A5_intron47/csv/`
  - Output: `plots_alphagenome_ISM_intron47/`
- `alphagenome_analysis_notebook_ISM_all_introns.ipynb`
  - Gene-agnostic postprocessing across **every intron of a gene**, consuming a
    full-gene batch ISM run (`../01_alphagenome_analysis/alphagenome_ISM_<GENE>/`,
    produced by `ISM/alphagenome_ISM_batch.ipynb`).
  - Config cell: set `GENE_NAME` (default `COL4A5`), `INPUT_ROOT`, `OUTPUT_DIR`,
    and optional `INTRON_FILTER`. Loads either many per-variant CSVs or a single
    aggregated CSV per intron (skips in-progress `.csv.tmp` files).
  - Applies the same pathogenicity cutoffs and spatial hotspot tests as the
    intron-6 workflow to each `intron_*` folder, then builds a cohort-level view
    comparing hotspot signal across introns: per-intron H/M burden (full and
    splice-adjacent-trimmed interior), deep-intronic hotspot **concentration**,
    and **focused ISM heatmaps** for significant (BH q < FDR) hotspot regions.
  - **Hotspot concentration (section 8.4):** default `CONCENTRATION_MODE = 'top_hotspot'`
    uses one representative interior hotspot per intron (the interval with the
    largest `n_HM_in_hotspot`), so concentration =
    `n_HM_in_representative_hotspot / n_HM_interior`. Set
    `CONCENTRATION_MODE = 'sum_hotspots'` to sum H/M across all significant
    interior hotspots instead. Outputs use suffix `_tophotspot` when in top-hotspot
    mode (e.g. `cohort_hotspot_concentration_tophotspot.png/svg/csv`).
  - Input: `../01_alphagenome_analysis/alphagenome_ISM_<GENE>/` (per-intron folders).
  - Output: `plots_alphagenome_ISM_<GENE>_all_introns/` (per-intron subfolders +
    `cohort_intron_summary.csv`, `all_introns_hotspot_intervals.csv`,
    `cohort_overview_HM_burden*.png/svg`, `cohort_hotspot_concentration*.png/svg`,
    and per-intron `hotspot_heatmaps/`).
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
  - Independent in-R reimplementation of motif MSEA, complementary to the Python/CLI preranked MSEA run.
  - Parameters cell: set `input_dir` and `output_dir` for the per-intron folder
    (e.g. `input_dir = plots_alphagenome_ISM_Intron6/` and
    `output_dir = plots_alphagenome_ISM_Intron6/complementary_motif_analysis/`).
  - Reads `motif_scores.rnk` and `motif_sets.gmt` from `input_dir`.
  - Uses `clusterProfiler::GSEA` (via `fgsea`), with `enrichplot`, `ggseqlogo`, and `pheatmap` for visualization.
  - Writes all outputs to `complementary_motif_analysis/` under the chosen intron folder:
    - `gsea_results_r.csv`, `gsea_significant_sets_r.csv`
    - `gsea_top_enriched_sets_r.csv`, `gsea_top_depleted_sets_r.csv`
    - `gsea_top_sets_nes_r.png`
    - `gsea_curve_enriched_<set>_r.png` and `gsea_curve_depleted_<set>_r.png` (top-N curves, default N = 5)
    - `leading_edge_membership_heatmap_r.png` and a leading-edge sequence-logo PNG
  - Tunable parameters in the first cell: `min_size`, `max_size`, `top_curves`, RNG seed.
  - Run once per hotspot intron (1, 4, 6, 30, 44, 49) before the cohort MSEA summary script.

## Scripts in this folder

- `motif_enrichment_pipeline.py`
  - Builds 7-mer motif scores and MSEA input files for a **single hotspot intron** per run.
  - Edit the PARAMETERS block: `INPUT_CSV`, `REFERENCE_SEQ`, `region_start`, and
    `OUTPUT_DIR`. Comment blocks document all six hotspot introns (1, 4, 6, 30, 44, 49)
    with genomic coordinates and reference sequences.
  - Input: `ism_variant_summary.csv` for the chosen intron (from
    `plots_alphagenome_ISM_COL4A5_all_introns/intron_<N>/`; for intron 6 motif/MSEA,
    prefer the all-introns table or the dedicated `plots_alphagenome_ISM_Intron6/` copy).
  - Builds a 7-mer (`Context7mer`) variant context using the hard-coded `REFERENCE_SEQ`
    and `region_start`; window size is controlled by `WINDOW` (5–8).
  - Ranks each 7-mer by the **mean of per-variant `Splice Site Max Quantile`**
    among variants that create that context. Configure source/label columns via
    `VARIANT_SCORE_COLUMN` and `AGG_SCORE_COLUMN`.
  - Outputs (written to `OUTPUT_DIR`, default `plots_alphagenome_ISM_Intron<N>/`):
    - `motif_scores_table.csv` — includes a `Creating Variants` column listing the
      unique `Variant ID`s that produced each motif (preserves motif-to-variant traceability).
    - `motif_scores.rnk` — two-column preranked MSEA input (`motif`, score).
    - `motif_sets.gmt` — motif sets built from non-positional 3-mers, simple
      homopolymer repeats (`X-repeat`), and `CG-core`/`AG-core` families.
    - `motif_sequences.fasta` — one record per unique motif, used by MEME / `ggseqlogo`.
  - These files feed both the CLI preranked MSEA workflow and the R notebook above.

- `plot_gsea_cohort_summary.py`
  - Aggregates per-intron R MSEA results into one cohort heatmap.
  - Reads `plots_alphagenome_ISM_Intron<N>/complementary_motif_analysis/gsea_results_r.csv`
    for each hotspot intron (default: 1, 4, 6, 30, 44, 49).
  - Rebuilds `motif_enrichment_cohort_summary/gsea_results_all_introns.csv`.
  - Plots **strict** significance (`p.adjust < 0.05`) in full red/blue; for introns
    **1 and 44 only**, nominally enriched motifs (`p.adjust < 0.25`) that fail the
    strict cutoff are shown in light gray.
  - Outputs:
    - `motif_enrichment_cohort_summary/gsea_significant_motifs_cohort.png`
    - `motif_enrichment_cohort_summary/gsea_significant_motifs_cohort.csv`
      (includes `display_status`: `significant` vs `nominal_gray`)
  - Run from this folder:
    ```bash
    python plot_gsea_cohort_summary.py
    python plot_gsea_cohort_summary.py --introns 1 4 6 30 44 49
    ```

- `plot_hotspot_motif_disruptions.py`
  - Counts **strict** motif disruptions among High/Moderate ISM variants in each
    hotspot window, using `MOTIF_RULES` from `ISM_motif_annotation.ipynb`
    (first-match-wins regex dictionary).
  - A variant counts as disrupting a motif only when the wild-type 7-mer matches a
    known rule **and** the mutant 7-mer no longer matches that same rule.
  - Hotspot reference windows are defined in `INTRON_CONFIGS` for introns 1, 4, 6, 30, 44, 49.
  - Also writes FASTA of reference and mutant 7-mer contexts for **all** High/Moderate
    variants (not only strict disruptions). Headers use `Variant ID`
    (e.g. `>chrX:108570691:A>G`); the following line is the 7-mer sequence.
  - Outputs under `motif_disruption_cohort_summary/`:
    - `hm_variant_motif_disruptions.csv` — per-variant disruption table
    - `hm_ref_7mer_motifs.fasta` — wild-type 7-mers for all H/M variants (all hotspots)
    - `hm_mut_7mer_motifs.fasta` — mutant 7-mers for all H/M variants (all hotspots)
    - `intron_<N>/hm_ref_7mer_motifs.fasta`, `intron_<N>/hm_mut_7mer_motifs.fasta` —
      same FASTAs split per hotspot intron
    - `motif_disruption_by_class.csv`, `motif_disruption_by_factor.csv`, `motif_disruption_by_impact.csv`
    - `motif_disruption_audit.txt` — strict vs lenient counts and caveats
    - `motif_disruption_by_class.png`, `motif_disruption_ISS_vs_other.png`,
      `motif_disruption_factor_heatmap.png`, `motif_disruption_ISS_by_impact.png`,
      `motif_disruption_fraction_known.png`
  - Run from this folder:
    ```bash
    python plot_hotspot_motif_disruptions.py
    python plot_hotspot_motif_disruptions.py --introns 1 6 44
    ```

## Input dependencies

Expected upstream inputs come from `../01_alphagenome_analysis/`:

- ISM input folders:
  - Intron 6: `../01_alphagenome_analysis/alphagenome_ISM_COL4A5_intron6/csv/`
  - Intron 47: `../01_alphagenome_analysis/alphagenome_ISM_COL4A5_intron47/csv/`
  - All introns (full-gene batch run): `../01_alphagenome_analysis/alphagenome_ISM_<GENE>/` (per-intron `intron_*/csv/` folders)
- Cohort input folder:
  - `../01_alphagenome_analysis/alphagenome_intron6_cohort/csv/`
- gnomAD input folder:
  - `../01_alphagenome_analysis/alphagenome_intron6_genomAD/csv/`

## Output folders and current artifacts

Output folder naming:
- Hotspot introns (1, 4, 6, 30, 44, 49): `plots_alphagenome_ISM_Intron<N>/` (capital `Intron`).
- Control intron 47: `plots_alphagenome_ISM_intron47/` (lowercase; separate control workflow).
- Legacy intron 6 single-window postprocessing: `plots_alphagenome_ISM/`.

- `plots_alphagenome_ISM/`
  - Legacy intron 6 postprocessing from `alphagenome_analysis_notebook_ISM.ipynb`:
    `ism_variant_summary.csv`, heatmaps, and motif files for the 160 nt single-window run.
- `plots_alphagenome_ISM_Intron<N>/` (hotspot introns 1, 4, 6, 30, 44, 49)
  - Per-intron motif pipeline outputs from `motif_enrichment_pipeline.py`:
    `motif_scores_table.csv`, `motif_scores.rnk`, `motif_sets.gmt`, `motif_sequences.fasta`
  - `complementary_motif_analysis/` — outputs of `run_complementary_motif_analysis.ipynb`
    (`gsea_results_r.csv`, `gsea_significant_sets_r.csv`, enrichment curves, etc.)
  - Intron 6 cohort scripts read and write `plots_alphagenome_ISM_Intron6/` (not the legacy `plots_alphagenome_ISM/` folder).
- `plots_alphagenome_ISM_<GENE>_all_introns/` (e.g. `plots_alphagenome_ISM_COL4A5_all_introns/`)
  - `cohort_intron_summary.csv`, `all_introns_hotspot_intervals.csv`
  - `cohort_overview_HM_burden.png/svg`, `cohort_overview_HM_burden_interior.png/svg`
  - `cohort_hotspot_concentration_tophotspot.png/svg/csv` (default concentration mode)
  - one `intron_*/` subfolder per intron with per-intron heatmaps, hotspot
    tables (`ism_variant_summary.csv`, `spatial_hotspots_HM_variants.csv`), and
    focused `hotspot_heatmaps/`
- `motif_enrichment_cohort_summary/`
  - `gsea_results_all_introns.csv`
  - `gsea_significant_motifs_cohort.png`, `gsea_significant_motifs_cohort.csv`
- `motif_disruption_cohort_summary/`
  - Per-variant and summary CSVs plus cohort disruption bar plots and heatmaps (see script section).
  - `hm_ref_7mer_motifs.fasta` / `hm_mut_7mer_motifs.fasta` — ref and mutant 7-mers
    for all High/Moderate variants across hotspots (`>chrX:pos:REF>ALT` headers).
  - `intron_<N>/` — per-hotspot copies of the same FASTAs.
- `plots_patients_cohort/`
  - `variant_scores_heatmap.svg`
  - `variant_genomic_positions.svg`
  - `variant_pathogenicity_summary.csv`
- `plots_alphagenome_genomAD/`
  - `variant_scores_heatmap.svg`
  - `variant_genomic_positions.svg`
  - `variant_pathogenicity_summary.csv`
- `plots_alphagenome_ISM_intron47/`
  - Intron 47 postprocessing outputs from `alphagenome_analysis_notebook_ISM_intron47.ipynb`.

## Motif enrichment (MSEA) output

- `motif_analysis.GseaPreranked.<timestamp>/`
  - Contains preranked MSEA report pages, `edb/` files, and enrichment tables generated from the `.rnk` / `.gmt` inputs.

## Recommended run order

1. Run upstream analyses in `../01_alphagenome_analysis/`.
2. Run `alphagenome_analysis_notebook_ISM.ipynb` (intron 6) and/or `alphagenome_analysis_notebook_ISM_intron47.ipynb` (intron 47). For a full-gene batch run, run `alphagenome_analysis_notebook_ISM_all_introns.ipynb` to postprocess every intron and build the cross-intron hotspot comparison.
3. Run `alphagenome_analysis_notebook_cohort.ipynb`.
4. Run `alphagenome_analysis_notebook_genomAD.ipynb`.
5. For each hotspot intron (1, 4, 6, 30, 44, 49):
   - Edit and run `motif_enrichment_pipeline.py` (set `INPUT_CSV`, `REFERENCE_SEQ`, `region_start`, `OUTPUT_DIR`).
   - Run `run_complementary_motif_analysis.ipynb` with matching `input_dir` / `output_dir`.
6. Run `plot_gsea_cohort_summary.py` to build the cohort MSEA NES heatmap.
7. Run `plot_hotspot_motif_disruptions.py` for strict ISS/ESE disruption counts and
   H/M ref/mutant 7-mer FASTAs across hotspots.
8. (Optional) Run `ISM_motif_annotation.ipynb` on a per-intron `motif_scores_table.csv`.
9. (Optional) run preranked MSEA using the generated `.rnk` and `.gmt`, then review outputs in `motif_analysis.GseaPreranked.*`.

## Environment

- Python notebooks/scripts: use the Conda environment defined in `coda/alphagenome-env.yml`. Notebook kernel: `alphagenome-env` (Python 3.11).
- `motif_enrichment_pipeline.py` additionally requires `biopython` (`from Bio.Seq import Seq`).
- `plot_gsea_cohort_summary.py` and `plot_hotspot_motif_disruptions.py` require `pandas`, `matplotlib`, and `numpy`.
- `run_complementary_motif_analysis.ipynb` runs on a separate R kernel `R (r-motif-gsea-env)` with Bioconductor packages `clusterProfiler`, `enrichplot`, `fgsea`, and CRAN packages `ggplot2`, `ggseqlogo`, `pheatmap`, `gridExtra` (install commands are commented in the first setup cell).
