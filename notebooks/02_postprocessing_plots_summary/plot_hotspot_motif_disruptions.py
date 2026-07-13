#!/usr/bin/env python3
"""
Count ISS and other known regulatory motifs disrupted by H/M ISM variants per hotspot.

Uses MOTIF_RULES from ISM_motif_annotation.ipynb (first-match-wins regex dictionary).

A variant counts as disrupting a motif only when:
  1. The wild-type 7-mer matches a known rule, AND
  2. The mutant 7-mer no longer matches that same rule (strict loss of motif).

Outputs: motif_disruption_cohort_summary/
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / "motif_disruption_cohort_summary"

HM_IMPACT = ("High Impact", "Moderate Impact")
WINDOW = 7
FLANK = WINDOW // 2


@dataclass(frozen=True)
class IntronConfig:
    intron: int
    region_start: int
    reference_seq: str
    input_csv: Path


INTRON_CONFIGS: dict[int, IntronConfig] = {
    1: IntronConfig(
        intron=1,
        region_start=108535962,
        reference_seq=(
            "CTGAATCAAAATTACTTAGAGATGTTTATTAAAATTTCAGATTCTGCTATTCAACTTCAGCTTACTGGAAGAATCTTTACGATGGGTTGCTGGTATGTAAATTTTTGTGAGGTTCTTTAGGTGACCCTAATTCA"
        ),
        input_csv=SCRIPT_DIR
        / "plots_alphagenome_ISM_COL4A5_all_introns/intron_1/ism_variant_summary.csv",
    ),
    4: IntronConfig(
        intron=4,
        region_start=108568001,
        reference_seq=(
            "TTTGTTTCTCATATTTTGTTTCATAATTTCTAGGGTATCTTCAGGACAGATTCTTAGAAGTGAAATGCTAGGTCAAAGGGCAAATGCATGGGTAATTTTGCTAGATATTTTCAAATTTCTCTCCATA"
        ),
        input_csv=SCRIPT_DIR
        / "plots_alphagenome_ISM_COL4A5_all_introns/intron_4/ism_variant_summary.csv",
    ),
    6: IntronConfig(
        intron=6,
        region_start=108570633,
        reference_seq=(
            "TAAACTTGATGTCTAGGCCACTTCCTTTCTCTCGGGACCTACTTTTTCCATGTGTAACAAGGTGGAGAGAAGGGTATTGGACTCACAAAGACACACAACAGTAGTAATTTTATTCTTTCAAACCTTCTGATGAAGTTGTTTCTAGGATTACCGTGGCATA"
        ),
        input_csv=SCRIPT_DIR / "plots_alphagenome_ISM_Intron6/ism_variant_summary.csv",
    ),
    30: IntronConfig(
        intron=30,
        region_start=108615585,
        reference_seq=(
            "TCTTTTGTCAACACGTAACCTAATGGGTTTGCCTCTATAGATAGGAAATCATATGGTTCCTCTAAATTTGTGTGAAAGCATATTGAGTTTCTGGATTCTGAGCTGTCTGGTGATGTGAATTCTCGTTATGTTAATCTAGGTAAGTACAGTA"
        ),
        input_csv=SCRIPT_DIR
        / "plots_alphagenome_ISM_COL4A5_all_introns/intron_30/ism_variant_summary.csv",
    ),
    44: IntronConfig(
        intron=44,
        region_start=108680168,
        reference_seq=(
            "ATGTTTCTTTCAAATGGCCAATGGGCAGGGGGTGTGTTCAAGCCACAGTTTATCTTCGCAGTGGTTTGACTTGATTGTACAAAGGGAGATCCTGTGACACTTTTACAACATCTTTCCAAGAACACAAGGTATGTGACATTTCCCTTACCCAACCACG"
        ),
        input_csv=SCRIPT_DIR
        / "plots_alphagenome_ISM_COL4A5_all_introns/intron_44/ism_variant_summary.csv",
    ),
    49: IntronConfig(
        intron=49,
        region_start=108688753,
        reference_seq=(
            "CTGACTATCCCCTTTGCCTCCATAAAGACATTTAAGAAGATAGTCTGGGCTCCCAGTTTCTGGGATCTCTATGATGAGAAGCTTCCTGAGATCTCATTGTCTGGGCCATTCATTCTGTGGAAATGTTTCTGCAGCAGCAGAGACCCCGGAGGTCAGCTTTACTCCCTATTTTTCCAAGGCATATGGCCACACCAGTTCTGGCCCACCAGGGTAAGGATATTCACTTAC"
        ),
        input_csv=SCRIPT_DIR
        / "plots_alphagenome_ISM_COL4A5_all_introns/intron_49/ism_variant_summary.csv",
    ),
}


def get_context(
    position: int,
    allele: str,
    region_start: int,
    reference_seq: str,
    flank: int = FLANK,
) -> str | None:
    idx = position - region_start
    start = idx - flank
    end = idx + flank + 1
    if start < 0 or end > len(reference_seq):
        return None
    context = list(reference_seq[start:end])
    context[flank] = allele.upper()
    return "".join(context).upper()


def load_hm_contexts(cfg: IntronConfig) -> pd.DataFrame:
    region_end = cfg.region_start + len(cfg.reference_seq) - 1
    df = pd.read_csv(cfg.input_csv)
    df = df[
        (df["Position"] >= cfg.region_start) & (df["Position"] <= region_end)
    ].copy()
    df = df[df["Impact Classification"].isin(HM_IMPACT)].copy()
    df["Context7mer"] = df.apply(
        lambda r: get_context(
            int(r["Position"]),
            str(r["Alt"]),
            cfg.region_start,
            cfg.reference_seq,
        ),
        axis=1,
    )
    return df[df["Context7mer"].notna()].copy()

# Mirrors ISM_motif_annotation.ipynb — ordered most-specific to least-specific.
MOTIF_RULES: list[tuple[str, str, str]] = [
    (r"CTAGG|TTAGG|TAGG", "hnRNP A1/A2", "ISS"),
    (r"TTCTT|TCTTT|TCTT", "PTBP1", "ISS"),
    (r"CTTC", "PTBP1/hnRNP A1", "ISS"),
    (r"TCTCTC|TCTCGG|TCTCG", "PTBP1", "ISS"),
    (r"GGGGG|GGGG|GGG", "hnRNP H/F", "ISS"),
    (r"TTTTT|TTTTTT", "TIA1/hnRNP C", "ISS"),
    (r"TTTT", "TIA1/hnRNP C", "ISS"),
    (r"GAAGAA|AAGAAG", "SRSF1", "ESE"),
    (r"GGAGG", "SRSF1", "ESE"),
    (r"AGTAAG", "SRSF2/SC35", "ESE"),
    (r"AGGAC", "SRSF5/SRp40", "ESE"),
    (r"[CT]CA[CT]", "NOVA1/2", "ISE/ISS"),
    (r"TGCATG|GCATG", "RBFOX1/2", "ISE/ISS"),
    (r"CG", "CpG-associated", "ISS"),
]

CLASS_ORDER = ["ISS", "ESE", "ISE/ISS"]
CLASS_COLORS = {"ISS": "#1b9e77", "ESE": "#d95f02", "ISE/ISS": "#7570b3"}
FIGURE_DPI = 300


def match_rule(seq: str) -> tuple[str | None, str, str]:
    for pattern, factor, cls in MOTIF_RULES:
        if re.search(pattern, seq):
            return pattern, factor, cls
    return None, "Unknown", "Unknown"


def motif_disrupted(ref_7mer: str, mut_7mer: str) -> bool:
    pattern, _factor, cls = match_rule(ref_7mer)
    if pattern is None or cls == "Unknown":
        return False
    return re.search(pattern, mut_7mer) is None


def add_ref_context(df: pd.DataFrame, cfg: IntronConfig) -> pd.DataFrame:
    out = df.copy()
    out["Reference7mer"] = out.apply(
        lambda r: get_context(int(r["Position"]), str(r["Ref"]), cfg.region_start, cfg.reference_seq),
        axis=1,
    )
    out = out[out["Reference7mer"].notna()].copy()

    patterns, factors, classes, disrupted = [], [], [], []
    for ref, mut in zip(out["Reference7mer"], out["Context7mer"]):
        pat, fac, cls = match_rule(ref)
        patterns.append(pat or "")
        factors.append(fac)
        classes.append(cls)
        disrupted.append(motif_disrupted(ref, mut))

    out["Matched_Pattern"] = patterns
    out["Regulatory_Factor"] = factors
    out["Regulatory_Class"] = classes
    out["motif_disrupted"] = disrupted
    return out


def build_variant_table(introns: list[int]) -> pd.DataFrame:
    parts = []
    for intron in introns:
        cfg = INTRON_CONFIGS[intron]
        hm = add_ref_context(load_hm_contexts(cfg), cfg)
        hm["intron"] = intron
        hm["impact"] = hm["Impact Classification"]
        parts.append(hm)
    return pd.concat(parts, ignore_index=True)


def disrupted_only(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["motif_disrupted"]].copy()


def summarize_by_intron(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    d = disrupted_only(df)
    for intron in sorted(df["intron"].unique()):
        n_total = len(df[df["intron"] == intron])
        sub = d[d["intron"] == intron]
        for cls in CLASS_ORDER:
            n_cls = int((sub["Regulatory_Class"] == cls).sum())
            rows.append({
                "intron": intron,
                "regulatory_class": cls,
                "n_variants_disrupting": n_cls,
                "fraction_of_HM": n_cls / n_total if n_total else 0.0,
            })
        rows.append({
            "intron": intron,
            "regulatory_class": "any_known",
            "n_variants_disrupting": len(sub),
            "fraction_of_HM": len(sub) / n_total if n_total else 0.0,
        })
    return pd.DataFrame(rows)


def summarize_by_factor(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for intron, sub in disrupted_only(df).groupby("intron"):
        for factor, count in sub["Regulatory_Factor"].value_counts().items():
            rows.append({
                "intron": intron,
                "regulatory_factor": factor,
                "n_variants_disrupting": int(count),
            })
    return pd.DataFrame(rows)


def summarize_by_impact(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    d = disrupted_only(df)
    for (intron, impact), sub in d.groupby(["intron", "impact"]):
        n_total = len(df[(df["intron"] == intron) & (df["impact"] == impact)])
        for cls in CLASS_ORDER + ["any_known"]:
            n_cls = len(sub) if cls == "any_known" else int((sub["Regulatory_Class"] == cls).sum())
            rows.append({
                "intron": intron,
                "impact": impact,
                "regulatory_class": cls,
                "n_variants_disrupting": n_cls,
                "n_HM_total": n_total,
            })
    return pd.DataFrame(rows)


def plot_class_counts(class_summary: pd.DataFrame, df: pd.DataFrame, path: Path) -> None:
    introns = sorted(class_summary["intron"].unique())
    x = np.arange(len(introns))
    width = 0.62
    bottom = np.zeros(len(introns))

    fig, ax = plt.subplots(figsize=(5.8, 3.2))
    for cls in CLASS_ORDER:
        vals = [
            int(class_summary.loc[
                (class_summary["intron"] == i) & (class_summary["regulatory_class"] == cls),
                "n_variants_disrupting",
            ].iloc[0])
            for i in introns
        ]
        ax.bar(x, vals, width, bottom=bottom, label=cls, color=CLASS_COLORS[cls],
               edgecolor="white", linewidth=0.6)
        for xi, val, base in zip(x, vals, bottom):
            if val >= 2:
                ax.text(xi, base + val / 2, str(val), ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
        bottom += np.array(vals)

    totals = [len(df[df["intron"] == i]) for i in introns]
    ymax = max(bottom) if len(bottom) else 1
    label_offset = max(1.8, ymax * 0.06)
    for xi, total, top in zip(x, totals, bottom):
        ax.text(xi, top + label_offset, f"n={total}", ha="center", va="bottom", fontsize=7, color="#444444")

    ax.set_xticks(x)
    ax.set_xticklabels([f"Intron {i}" for i in introns], fontsize=9)
    ax.set_ylabel("H/M spliceogenic variants\ncausing motif loss", fontsize=9)
    ax.set_title(
        "Strict motif disruptions per hotspot\n(wt 7-mer matched rule; mut 7-mer does not)",
        fontsize=9,
    )
    ax.set_ylim(0, ymax + label_offset * 2.2)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="Regulatory class", frameon=False, loc="upper right", fontsize=8, title_fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_iss_vs_other(class_summary: pd.DataFrame, path: Path) -> None:
    introns = sorted(class_summary["intron"].unique())
    iss, other = [], []
    for i in introns:
        sub = class_summary[(class_summary["intron"] == i) & (class_summary["regulatory_class"].isin(CLASS_ORDER))]
        iss.append(int(sub.loc[sub["regulatory_class"] == "ISS", "n_variants_disrupting"].sum()))
        other.append(int(sub.loc[sub["regulatory_class"] != "ISS", "n_variants_disrupting"].sum()))
    x = np.arange(len(introns))
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ax.bar(x, iss, 0.6, label="ISS", color=CLASS_COLORS["ISS"])
    ax.bar(x, other, 0.6, bottom=iss, label="ESE + ISE/ISS", color="#fdae6b")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Intron {i}" for i in introns])
    ax.set_ylabel("Variants with strict motif loss")
    ax.set_title("ISS vs other motif disruptions (strict)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_factor_heatmap(factor_summary: pd.DataFrame, path: Path) -> None:
    if factor_summary.empty:
        return
    factors = sorted(factor_summary["regulatory_factor"].unique())
    introns = sorted(factor_summary["intron"].unique())
    mat = np.zeros((len(factors), len(introns)))
    for j, intron in enumerate(introns):
        for i, factor in enumerate(factors):
            row = factor_summary[
                (factor_summary["intron"] == intron) & (factor_summary["regulatory_factor"] == factor)
            ]
            mat[i, j] = int(row["n_variants_disrupting"].iloc[0]) if len(row) else 0

    n_rows, n_cols = len(factors), len(introns)
    fig, ax = plt.subplots(figsize=(5.6, 0.28 * n_rows + 1.6))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([f"Intron {i}" for i in introns], fontsize=9)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(factors, fontsize=8)
    vmax = mat.max() if mat.size else 1
    for i in range(n_rows):
        for j in range(n_cols):
            val = int(mat[i, j])
            if val > 0:
                text_color = "white" if val > vmax * 0.45 else "#333333"
                ax.text(j, i, str(val), ha="center", va="center", fontsize=8,
                        color=text_color, fontweight="bold")
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("Strictly disrupted RBP / motif per hotspot", fontsize=9, pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("# variants", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_by_impact(impact_summary: pd.DataFrame, path: Path) -> None:
    introns = sorted(impact_summary["intron"].unique())
    impacts = ["High Impact", "Moderate Impact"]
    x = np.arange(len(introns))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for k, impact in enumerate(impacts):
        vals = []
        for intron in introns:
            row = impact_summary[
                (impact_summary["intron"] == intron)
                & (impact_summary["impact"] == impact)
                & (impact_summary["regulatory_class"] == "ISS")
            ]
            vals.append(int(row["n_variants_disrupting"].iloc[0]) if len(row) else 0)
        color = "#d62728" if impact == "High Impact" else "#ff7f0e"
        ax.bar(x + (k - 0.5) * width, vals, width, label=impact.replace(" Impact", ""), color=color, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Intron {i}" for i in introns])
    ax.set_ylabel("ISS motif loss (strict)")
    ax.set_title("ISS disruptions — High vs Moderate")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_fraction_known(class_summary: pd.DataFrame, path: Path) -> None:
    introns = sorted(class_summary["intron"].unique())
    fracs = [
        float(class_summary.loc[
            (class_summary["intron"] == i) & (class_summary["regulatory_class"] == "any_known"),
            "fraction_of_HM",
        ].iloc[0])
        for i in introns
    ]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar([f"Intron {i}" for i in introns], [f * 100 for f in fracs], color="#4c72b0", alpha=0.85)
    ax.set_ylabel("% H/M with strict motif loss")
    ax.set_ylim(0, max(50, max(fracs) * 100 + 10))
    ax.set_title("Fraction of H/M variants with strict motif disruption")
    for bar, frac in zip(bars, fracs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                f"{frac * 100:.0f}%", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_audit_notes(df: pd.DataFrame, path: Path) -> None:
    lenient = (df["Regulatory_Class"] != "Unknown").sum()
    strict = df["motif_disrupted"].sum()
    ggg = int((df["Matched_Pattern"] == r"GGGGG|GGGG|GGG").sum())
    ggg_strict = int(((df["Matched_Pattern"] == r"GGGGG|GGGG|GGG") & df["motif_disrupted"]).sum())

    lines = [
        "Motif disruption audit notes",
        "=" * 50,
        "",
        "Annotation method: ISM_motif_annotation.ipynb MOTIF_RULES (regex, first match).",
        "",
        "STRICT disruption criterion:",
        "  wt 7-mer matches a rule AND mutant 7-mer no longer matches that same rule.",
        "",
        f"Lenient count (wt has any known motif): {lenient}",
        f"Strict count (motif actually lost in mut): {strict}",
        "",
        "Caveats:",
        "  • GGG/GGGG rule (hnRNP H/F) is very permissive — any 3+ G run in 7-mer.",
        f"    Variants with wt G-run match: {ggg}; strict G-run loss: {ggg_strict}",
        "  • CG rule only applies when no higher-priority rule matched first.",
        "  • Regex dictionary is heuristic — not experimentally validated per variant.",
        "  • Variants with no wt motif match are NOT counted as disruptions.",
        "",
        "Pattern breakdown (wt matches, all H/M):",
    ]
    for pat, n in df["Matched_Pattern"].value_counts().items():
        if pat:
            n_strict = int(((df["Matched_Pattern"] == pat) & df["motif_disrupted"]).sum())
            lines.append(f"  {pat}: wt={n}, strict_loss={n_strict}")
    path.write_text("\n".join(lines) + "\n")


def process(introns: list[int]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df = build_variant_table(introns)
    df.to_csv(OUTPUT_DIR / "hm_variant_motif_disruptions.csv", index=False)

    class_summary = summarize_by_intron(df)
    factor_summary = summarize_by_factor(df)
    impact_summary = summarize_by_impact(df)
    class_summary.to_csv(OUTPUT_DIR / "motif_disruption_by_class.csv", index=False)
    factor_summary.to_csv(OUTPUT_DIR / "motif_disruption_by_factor.csv", index=False)
    impact_summary.to_csv(OUTPUT_DIR / "motif_disruption_by_impact.csv", index=False)
    write_audit_notes(df, OUTPUT_DIR / "motif_disruption_audit.txt")

    plot_class_counts(class_summary, df, OUTPUT_DIR / "motif_disruption_by_class.png")
    plot_iss_vs_other(class_summary, OUTPUT_DIR / "motif_disruption_ISS_vs_other.png")
    plot_factor_heatmap(factor_summary, OUTPUT_DIR / "motif_disruption_factor_heatmap.png")
    plot_by_impact(impact_summary, OUTPUT_DIR / "motif_disruption_ISS_by_impact.png")
    plot_fraction_known(class_summary, OUTPUT_DIR / "motif_disruption_fraction_known.png")

    print(f"Wrote outputs to {OUTPUT_DIR}\n")
    print("STRICT disruptions per hotspot (ISS / ESE / ISE-ISS / any):")
    for intron in sorted(introns):
        sub = class_summary[class_summary["intron"] == intron]
        iss = int(sub.loc[sub["regulatory_class"] == "ISS", "n_variants_disrupting"].iloc[0])
        ese = int(sub.loc[sub["regulatory_class"] == "ESE", "n_variants_disrupting"].iloc[0])
        ise = int(sub.loc[sub["regulatory_class"] == "ISE/ISS", "n_variants_disrupting"].iloc[0])
        any_k = int(sub.loc[sub["regulatory_class"] == "any_known", "n_variants_disrupting"].iloc[0])
        total = len(df[df["intron"] == intron])
        print(f"  Intron {intron}: ISS={iss}, ESE={ese}, ISE/ISS={ise}, any={any_k}/{total} ({any_k/total*100:.0f}%)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot strict motif disruptions per hotspot.")
    parser.add_argument("--introns", type=int, nargs="+", default=sorted(INTRON_CONFIGS), choices=sorted(INTRON_CONFIGS))
    return parser.parse_args()


def main() -> None:
    process(parse_args().introns)


if __name__ == "__main__":
    main()
