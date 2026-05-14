#!/usr/bin/env python3
"""
Pipeline to identify splicing motifs from in‑silico mutagenesis data
and prepare input for GSEA + MEME motif visualization.

Author: Hassan Saei
Date: 2026-03-29
"""

import os

import pandas as pd
from Bio.Seq import Seq

# ===============================
# PARAMETERS
# ===============================
INPUT_CSV = "plots_alphagenome_ISM/ism_variant_summary.csv"

# genomic sequence defining our tested region
REFERENCE_SEQ = (
    "TAAACTTGATGTCTAGGCCACTTCCTTTCTCTCGGGACCTACTTTTTCCATGTGTAACAAGGTGGAGAGAAGGGTATTGGACTCACAAAGACACACAACAGTAGTAATTTTATTCTTTCAAACCTTCTGATGAAGTTGTTTCTAGGATTACCGTGGCATA"
)
region_start = 108570633

OUTPUT_DIR = "plots_alphagenome_ISM"
OUTPUT_TABLE = os.path.join(OUTPUT_DIR, "motif_scores_table.csv")
OUTPUT_RNK = os.path.join(OUTPUT_DIR, "motif_scores.rnk")
OUTPUT_GMT = os.path.join(OUTPUT_DIR, "motif_sets.gmt")
OUTPUT_FASTA = os.path.join(OUTPUT_DIR, "motif_sequences.fasta")
WINDOW = 7  # size for motif scoring (can be 5–8)

# Per-variant column in ism_variant_summary.csv (cohort-aligned: splice-site usage)
VARIANT_SCORE_COLUMN = "Splice Site Max Quantile"
# After groupby(Context7mer), we take the mean of that column per motif:
AGG_SCORE_COLUMN = "Splice Site Max Quantile"

# ===============================
# HELPER FUNCTION
# ===============================
def get_context(position, ref, alt, flank=3):
    """
    Extract (2*flank+1)-nt context around a mutation within the reference sequence
    and substitute the ALT base in place of REF.
    """
    idx = position - region_start          # 0‑based coordinate relative to seq start
    start = max(0, idx - flank)
    end   = min(len(REFERENCE_SEQ), idx + flank + 1)
    context = list(REFERENCE_SEQ[start:end])
    if 0 <= idx - start < len(context):
        context[idx - start] = alt.upper()
    return "".join(context).upper()


def join_unique_variant_ids(values):
    """Return unique variant IDs in appearance order."""
    unique_ids = dict.fromkeys(v for v in values if pd.notna(v))
    return ";".join(str(v) for v in unique_ids)

# ===============================
# LOAD VARIANTS
# ===============================

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(INPUT_CSV)
if VARIANT_SCORE_COLUMN not in df.columns:
    raise KeyError(
        f"Missing column {VARIANT_SCORE_COLUMN!r} in {INPUT_CSV}. "
        "Regenerate ism_variant_summary.csv from alphagenome_analysis_notebook_ISM.ipynb."
    )

# optional: keep all impact levels so ranking covers full spectrum
# df = df[df["Impact Classification"].isin(["High Impact",
#                                           "Moderate Impact",
#                                           "Low Impact"])]

# context around each variant (7‑mer window gives flexibility)
df["Context7mer"] = df.apply(
    lambda r: get_context(r["Position"], r["Ref"], r["Alt"], flank=3),
    axis=1,
)
# shorter fixed window for standard motif analysis (e.g. 7‑mer)
df["Context7mer"] = df.apply(
    lambda r: get_context(r["Position"], r["Ref"], r["Alt"], flank=WINDOW // 2),
    axis=1,
)

# ===============================
# COMPUTE & RANK MOTIF SCORES
# ===============================

motif_scores = (
    df.groupby("Context7mer")
      .agg(
          **{
              AGG_SCORE_COLUMN: (VARIANT_SCORE_COLUMN, "mean"),
              "Creating Variants": ("Variant ID", join_unique_variant_ids),
          }
      )
      .sort_values(by=AGG_SCORE_COLUMN, ascending=False)
      .reset_index()
)
motif_scores.to_csv(OUTPUT_TABLE, index=False)

# create GSEA .rnk file  (motif\tScore)
motif_scores[["Context7mer", AGG_SCORE_COLUMN]].to_csv(
    OUTPUT_RNK, sep="\t", index=False, header=False
)

# ===============================
# BUILD MOTIF SETS (.gmt)
# ===============================

def make_motif_sets():
    motifs = motif_scores["Context7mer"].tolist()
    motif_sets = {}

    for m in motifs:
        # non‑positional 3‑mers
        for i in range(len(m) - 2):
            triplet = m[i:i+3]
            motif_sets.setdefault(triplet, []).append(m)
        # simple repeats (AAAAA, CCCCC, etc.)
        if len(set(m)) == 1:
            motif_sets.setdefault(f"{m[0]}-repeat", []).append(m)
        # core‑dinucleotide families
        if "CG" in m:
            motif_sets.setdefault("CG-core", []).append(m)
        if "AG" in m:
            motif_sets.setdefault("AG-core", []).append(m)

    with open(OUTPUT_GMT, "w") as g:
        for name, members in motif_sets.items():
            members_unique = sorted(set(members))
            g.write(f"{name}\tNA\t" + "\t".join(members_unique) + "\n")

make_motif_sets()

# ===============================
# FASTA FOR MEME / ggseqlogo
# ===============================

with open(OUTPUT_FASTA, "w") as f:
    for _, row in motif_scores.iterrows():
        motif = row["Context7mer"]
        f.write(f">{motif}\n{motif}\n")

print("✅ Files generated:")
print(f"  - Motif table   : {OUTPUT_TABLE}")
print(f"  - Motif ranking : {OUTPUT_RNK}")
print(f"  - Motif sets    : {OUTPUT_GMT}")
print(f"  - FASTA         : {OUTPUT_FASTA}")