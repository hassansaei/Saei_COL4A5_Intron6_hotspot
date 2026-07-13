#!/usr/bin/env python3
"""
Summarize GSEA significant motif sets across COL4A5 hotspot introns in one plot.

Strict significance: p.adjust < 0.05 (full color).
For introns 1 and 44, nominally enriched motifs (p.adjust < 0.25) are shown in
light gray when they fail the strict cutoff, so enrichment can be discussed
without claiming FDR significance.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "motif_enrichment_cohort_summary"
COHORT_CSV = OUTPUT_DIR / "gsea_results_all_introns.csv"

DEFAULT_INTRONS = [1, 4, 6, 30, 44, 49]
PADJ_STRICT = 0.05
PADJ_NOMINAL = 0.25
NOMINAL_DISPLAY_INTRONS = {1, 44}
TOP_ENRICHED = 3
TOP_DEPLETED = 1
FIGURE_DPI = 300
GRAY_FILL = "#d9d9d9"
GRAY_TEXT = "#555555"


def intron_dir(intron: int) -> Path:
    return SCRIPT_DIR / f"plots_alphagenome_ISM_Intron{intron}"


def load_per_intron_gsea(intron: int) -> pd.DataFrame:
    path = intron_dir(intron) / "complementary_motif_analysis/gsea_results_r.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing GSEA results for intron {intron}: {path}")
    df = pd.read_csv(path)
    df["intron"] = intron
    df["intron_label"] = f"Intron {intron}"
    return df


def load_per_intron_significant(intron: int) -> pd.DataFrame:
    path = intron_dir(intron) / "complementary_motif_analysis/gsea_significant_sets_r.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing significant GSEA sets for intron {intron}: {path}")
    return pd.read_csv(path)


def rebuild_cohort_csv(introns: list[int], path: Path = COHORT_CSV) -> pd.DataFrame:
    parts = [load_per_intron_gsea(i) for i in introns]
    df = pd.concat(parts, ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df


def select_summary_motifs(
    df: pd.DataFrame,
    introns: list[int],
    padj_strict: float,
    top_enriched: int = TOP_ENRICHED,
    top_depleted: int = TOP_DEPLETED,
) -> tuple[list[str], dict[str, set[int]]]:
    """
    Return motif row order and which introns nominated each motif.

    Rows include:
      - motifs with p.adjust < padj_strict in any intron
      - top enriched/depleted sets from each intron's nominal significant list
    """
    df = df[df["intron"].isin(introns)].copy()
    strict_ids = sorted(df.loc[df["p.adjust"] < padj_strict, "ID"].drop_duplicates().tolist())

    selected: list[str] = []
    sources: dict[str, set[int]] = {}
    seen: set[str] = set()

    def add_motifs(motifs: list[str], intron: int) -> None:
        for motif in motifs:
            sources.setdefault(motif, set()).add(intron)
            if motif not in seen:
                selected.append(motif)
                seen.add(motif)

    for intron in introns:
        sig = load_per_intron_significant(intron)
        enriched = sig[sig["NES"] > 0].nlargest(top_enriched, "NES")["ID"].tolist()
        depleted = sig[sig["NES"] < 0].nsmallest(top_depleted, "NES")["ID"].tolist()
        add_motifs(enriched + depleted, intron)

    for motif in strict_ids:
        if motif not in seen:
            selected.append(motif)
            seen.add(motif)
        sources.setdefault(motif, set())
        for intron in introns:
            row = df[(df["ID"] == motif) & (df["intron"] == intron)]
            if len(row) and float(row["p.adjust"].iloc[0]) < padj_strict:
                sources[motif].add(intron)

    if "CG-core" in seen:
        selected = ["CG-core"] + [m for m in selected if m != "CG-core"]
    return selected, sources


def build_matrices(
    df: pd.DataFrame,
    motifs: list[str],
    introns: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    nes = df.pivot(index="ID", columns="intron", values="NES").reindex(motifs)
    padj = df.pivot(index="ID", columns="intron", values="p.adjust").reindex(motifs)
    return nes[introns], padj[introns]


def cell_is_nominal_gray(
    intron: int,
    motif: str,
    padj_val: float,
    sources: dict[str, set[int]],
    padj_strict: float,
    padj_nominal: float,
) -> bool:
    if pd.isna(padj_val):
        return False
    if padj_val < padj_strict:
        return False
    if intron not in NOMINAL_DISPLAY_INTRONS:
        return False
    if padj_val >= padj_nominal:
        return False
    return intron in sources.get(motif, set())


def plot_gsea_cohort_summary(
    nes: pd.DataFrame,
    padj: pd.DataFrame,
    sources: dict[str, set[int]],
    path: Path,
    padj_strict: float = PADJ_STRICT,
    padj_nominal: float = PADJ_NOMINAL,
) -> None:
    motifs = list(nes.index)
    introns = list(nes.columns)
    data = nes.to_numpy(dtype=float)
    padj_arr = padj.to_numpy(dtype=float)

    strict = padj_arr < padj_strict
    gray = np.zeros_like(strict, dtype=bool)
    for i, motif in enumerate(motifs):
        for j, intron in enumerate(introns):
            gray[i, j] = cell_is_nominal_gray(
                intron, motif, padj_arr[i, j], sources, padj_strict, padj_nominal,
            )

    display = strict | gray
    masked_strict = np.ma.masked_where(~strict, data)

    n_rows, n_cols = len(motifs), len(introns)
    fig_w = 0.50 * n_cols + 2.2
    fig_h = 0.20 * n_rows + 1.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 0.08], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    right = gs[0, 1].subgridspec(2, 1, height_ratios=[1.0, 2.4], hspace=0.10)
    leg_ax = fig.add_subplot(right[0, 0])
    cax = fig.add_subplot(right[1, 0])

    vmax = max(1.8, float(np.nanmax(np.abs(data[strict])))) if strict.any() else 1.8
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="white")

    ax.set_facecolor("white")
    for i in range(n_rows):
        for j in range(n_cols):
            if gray[i, j]:
                ax.add_patch(mpatches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    facecolor=GRAY_FILL, edgecolor="white", linewidth=1.0,
                ))

    im = ax.imshow(masked_strict, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([f"Intron {i}" for i in introns], fontsize=9)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(motifs, fontsize=8)
    ax.set_title("GSEA motif per intron hotspot", fontsize=9, pad=8)

    for i in range(n_rows):
        for j in range(n_cols):
            if not display[i, j]:
                continue
            val = data[i, j]
            if strict[i, j]:
                text_color = "white" if abs(val) > vmax * 0.55 else "#222222"
            else:
                text_color = GRAY_TEXT
            ax.text(
                j, i, f"{val:.2f}", ha="center", va="center",
                fontsize=7, color=text_color, fontweight="bold",
            )

    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, cax=cax)
    cax.set_box_aspect(12)
    cbar.set_label("NES", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    legend_handles = [
        mpatches.Patch(facecolor="#b2182b", edgecolor="none", label=f"p.adj < {padj_strict:g} (enriched)"),
        mpatches.Patch(facecolor="#2166ac", edgecolor="none", label=f"p.adj < {padj_strict:g} (depleted)"),
        mpatches.Patch(
            facecolor=GRAY_FILL, edgecolor="#aaaaaa",
            label=f"Nominal (introns 1 & 44, p.adj ≥ {padj_strict:g})",
        ),
    ]
    leg_ax.axis("off")
    leg_ax.legend(handles=legend_handles, frameon=False, fontsize=7, loc="center left")

    fig.subplots_adjust(left=0.12, right=0.90, top=0.92, bottom=0.13)
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_significant_table(
    nes: pd.DataFrame,
    padj: pd.DataFrame,
    sources: dict[str, set[int]],
    path: Path,
    padj_strict: float,
    padj_nominal: float,
) -> None:
    rows = []
    for motif in nes.index:
        for intron in nes.columns:
            p = padj.loc[motif, intron]
            n = nes.loc[motif, intron]
            if pd.isna(p):
                continue
            if p < padj_strict:
                status = "significant"
            elif cell_is_nominal_gray(intron, motif, p, sources, padj_strict, padj_nominal):
                status = "nominal_gray"
            else:
                continue
            rows.append({
                "intron": intron,
                "motif_set": motif,
                "NES": n,
                "p.adjust": p,
                "direction": "enriched" if n > 0 else "depleted",
                "display_status": status,
            })
    pd.DataFrame(rows).sort_values(["intron", "p.adjust"]).to_csv(path, index=False)


def process(
    introns: list[int] | None = None,
    padj_strict: float = PADJ_STRICT,
    padj_nominal: float = PADJ_NOMINAL,
) -> None:
    introns = introns or DEFAULT_INTRONS
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = rebuild_cohort_csv(introns)
    motifs, sources = select_summary_motifs(df, introns, padj_strict)
    nes, padj = build_matrices(df, motifs, introns)

    plot_gsea_cohort_summary(
        nes, padj, sources,
        OUTPUT_DIR / "gsea_significant_motifs_cohort.png",
        padj_strict=padj_strict,
        padj_nominal=padj_nominal,
    )
    write_significant_table(
        nes, padj, sources,
        OUTPUT_DIR / "gsea_significant_motifs_cohort.csv",
        padj_strict=padj_strict,
        padj_nominal=padj_nominal,
    )

    strict_n = int((padj < padj_strict).sum().sum())
    gray_n = sum(
        cell_is_nominal_gray(int(intron), motif, padj.loc[motif, intron], sources, padj_strict, padj_nominal)
        for motif in nes.index
        for intron in nes.columns
    )
    print(f"Rebuilt {COHORT_CSV}")
    print(f"Wrote {OUTPUT_DIR / 'gsea_significant_motifs_cohort.png'}")
    print(f"Motif sets in plot: {len(motifs)} | strict: {strict_n} | gray (introns 1 & 44): {gray_n}")
    for intron in introns:
        n_strict = int((padj[intron] < padj_strict).sum())
        n_gray = sum(
            cell_is_nominal_gray(intron, motif, padj.loc[motif, intron], sources, padj_strict, padj_nominal)
            for motif in nes.index
        )
        print(f"  Intron {intron}: strict={n_strict}, gray={n_gray}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot cohort GSEA significant motifs.")
    parser.add_argument("--introns", type=int, nargs="+", default=DEFAULT_INTRONS)
    parser.add_argument("--padj-strict", type=float, default=PADJ_STRICT)
    parser.add_argument("--padj-nominal", type=float, default=PADJ_NOMINAL)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    process(
        introns=args.introns,
        padj_strict=args.padj_strict,
        padj_nominal=args.padj_nominal,
    )


if __name__ == "__main__":
    main()
