#!/usr/bin/env python3
"""
Pipeline to identify splicing motifs from in‑silico mutagenesis data
and prepare input for GSEA + MEME motif visualization.

Author: Hassan Saei
Date: 2026-03-29
"""

import os

import pandas as pd

# ===============================
# PARAMETERS
# ===============================
# To analyse a different intron, change REFERENCE_SEQ + region_start together
# (and point INPUT_CSV / OUTPUT_DIR at the matching intron). REFERENCE_SEQ MUST be
# the genomic plus-strand sequence whose FIRST base sits exactly at region_start,
# because variants are mapped with idx = Position - region_start. The reference
# alignment check below verifies this every run.
#
# Intron 6:  TAAACTTGATGTCTAGGCCACTTCCTTTCTCTCGGGACCTACTTTTTCCATGTGTAACAAGGTGGAGAGAAGGGTATTGGACTCACAAAGACACACAACAGTAGTAATTTTATTCTTTCAAACCTTCTGATGAAGTTGTTTCTAGGATTACCGTGGCATA
#            region_start = 108570633
# Intron 44: ATGTTTCTTTCAAATGGCCAATGGGCAGGGGGTGTGTTCAAGCCACAGTTTATCTTCGCAGTGGTTTGACTTGATTGTACAAAGGGAGATCCTGTGACACTTTTACAACATCTTTCCAAGAACACAAGGTATGTGACATTTCCCTTACCCAACCACG
#            region_start = 108680168
#
# Intron 30: TCTTTTGTCAACACGTAACCTAATGGGTTTGCCTCTATAGATAGGAAATCATATGGTTCCTCTAAATTTGTGTGAAAGCATATTGAGTTTCTGGATTCTGAGCTGTCTGGTGATGTGAATTCTCGTTATGTTAATCTAGGTAAGTACAGTA
#            region_start = 108615585
#
# Intron 49: CTGACTATCCCCTTTGCCTCCATAAAGACATTTAAGAAGATAGTCTGGGCTCCCAGTTTCTGGGATCTCTATGATGAGAAGCTTCCTGAGATCTCATTGTCTGGGCCATTCATTCTGTGGAAATGTTTCTGCAGCAGCAGAGACCCCGGAGGTCAGCTTTACTCCCTATTTTTCCAAGGCATATGGCCACACCAGTTCTGGCCCACCAGGGTAAGGATATTCACTTAC
#            region_start = 108688753
#
# Intron 4: TTTGTTTCTCATATTTTGTTTCATAATTTCTAGGGTATCTTCAGGACAGATTCTTAGAAGTGAAATGCTAGGTCAAAGGGCAAATGCATGGGTAATTTTGCTAGATATTTTCAAATTTCTCTCCATA
#            region_start = 108568001

INPUT_CSV = "plots_alphagenome_ISM_COL4A5_all_introns/intron_4/ism_variant_summary.csv"

REFERENCE_SEQ = (
    "TTTGTTTCTCATATTTTGTTTCATAATTTCTAGGGTATCTTCAGGACAGATTCTTAGAAGTGAAATGCTAGGTCAAAGGGCAAATGCATGGGTAATTTTGCTAGATATTTTCAAATTTCTCTCCATA"
)
region_start = 108568001

OUTPUT_DIR = "plots_alphagenome_ISM_Intron4"
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
def get_context(position, alt, flank=3):
    """
    Extract a (2*flank+1)-nt window centred on `position` from REFERENCE_SEQ and
    substitute the ALT base at the centre.

    Returns None when the full window is not contained within REFERENCE_SEQ, i.e.
    when the variant lies outside (or at the very edge of) the tested region. Such
    variants are dropped rather than silently yielding empty/truncated motifs.
    Note: idx assumes REFERENCE_SEQ[0] is the genomic plus-strand base at
    region_start (see check_reference_alignment).
    """
    idx = position - region_start          # 0-based coordinate relative to seq start
    start = idx - flank
    end = idx + flank + 1
    if start < 0 or end > len(REFERENCE_SEQ):
        return None
    context = list(REFERENCE_SEQ[start:end])
    context[flank] = alt.upper()           # centre of the window (idx - start == flank)
    return "".join(context).upper()


def check_reference_alignment(frame):
    """
    Verify REFERENCE_SEQ / region_start are consistent with the CSV: for every
    variant that falls inside the window, REFERENCE_SEQ[idx] must equal its Ref
    base. Raises if not, which catches a wrong sequence, start coordinate, or strand.
    """
    upper = region_start + len(REFERENCE_SEQ)
    inside = frame[(frame["Position"] >= region_start) & (frame["Position"] < upper)]
    mismatches = [
        (int(pos), REFERENCE_SEQ[int(pos) - region_start], ref)
        for pos, ref in zip(inside["Position"], inside["Ref"])
        if REFERENCE_SEQ[int(pos) - region_start] != ref
    ]
    n_inside = len(inside)
    print(
        f"Reference alignment: {n_inside - len(mismatches)}/{n_inside} "
        f"in-window Ref bases match REFERENCE_SEQ (region {region_start}-{upper - 1})"
    )
    if mismatches:
        raise ValueError(
            f"{len(mismatches)} in-window variants disagree with REFERENCE_SEQ "
            f"(e.g. (pos, seq_base, csv_ref)={mismatches[:5]}). "
            "Check REFERENCE_SEQ, region_start, and strand."
        )


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

# sanity check that REFERENCE_SEQ + region_start line up with this CSV
check_reference_alignment(df)

# WINDOW-mer context around each variant; None for variants outside the region
df["Context7mer"] = df.apply(
    lambda r: get_context(r["Position"], r["Alt"], flank=WINDOW // 2),
    axis=1,
)

# keep only variants whose full context fits inside the tested region
n_before = len(df)
df = df[df["Context7mer"].notna()].copy()
print(
    f"Kept {len(df)}/{n_before} variants within the "
    f"{len(REFERENCE_SEQ)}-nt tested region"
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