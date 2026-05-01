# Targeted RNA-seq — COL4A5 intron 6 sashimi analysis

This folder holds the targeted RNA-seq sample sheets, palette, and (in the
**gitignored** `bams/` subfolder) the per-replicate BAM files + `.bai`
indexes used to inspect splicing around the **COL4A5 intron 6 hotspot**
(cryptic / pseudoexon inclusion) on the human reference **hg19 / GRCh37**.

> **BAM files on Zenodo.** The 21 BAMs and their `.bai` indexes are not
> tracked in git; they are archived on Zenodo:
> [**10.5281/zenodo.19854065**](https://doi.org/10.5281/zenodo.19854065).
> Download the archive and place the files in `Targeted_rnaseq/bams/`
> (keep the original filenames) so the absolute paths in `samples.tsv`
> resolve before running ggsashimi.

### Folder layout

```
Targeted_rnaseq/
├── README.md           # this file
├── samples.tsv         # ggsashimi -b input: id<TAB>abs_bam_path<TAB>color_factor  (21 rows, committed)
├── palette.tsv         # ggsashimi -P input: one hex per line (4 rows, committed)
└── bams/               # gitignored — not pushed to GitHub
    ├── *.bam           # 21 BAM files
    └── *.bam.bai       # 21 index files
```

---

## 1. Sample design
Three biological groups, each with three replicates (R1–R3) where applicable.
Emetine (NMD inhibitor) was applied to allow detection of unstable
NMD-targeted transcripts.

| Group            | Treatment                          | Replicates             | Suffix  |
|------------------|------------------------------------|------------------------|---------|
| `AS_015`         | patient cells (Alport AS_015)      | R1, R2, R3 × En/Ep     | `_En` = emetine-negative, `_Ep` = emetine-positive |
| `AS_015_ASO`     | patient cells (AS_015) + ASO       | R1, R2, R3             | no En/Ep split |
| `control_1_56`   | control cell URECs line                  | R1, R2, R3 × En/Ep     | same    |
| `control_1_29`   | control cell URECs line                  | R1, R2, R3 × En/Ep     | same    |

Total: **21 BAM files** (+ `.bai` indexes).

### Rename history
Originally labelled as:

| Original prefix | Renamed to         | Reason |
|-----------------|--------------------|--------|
| `1_29_*`        | `AS_015_*`         | this cell line is the patient sample (Alport ID AS_015), not a control |
| `MW_*`          | `control_1_29_*`   | re-assigned to control group |
| `1_56_*`        | `control_1_56_*`   | explicit control labelling |
| `ASO_[1-3]`     | `AS_015_ASO_R[1-3]` | ASO samples are AS_015 patient cells treated with ASO; replicate suffix harmonised to R1/R2/R3 |

Applied to both `.bam` and `.bam.bai` files on 2026-04-21.
The `ASO_1_.bam` file (note trailing underscore) was also normalised to
`AS_015_ASO_1.bam` during the rename.

---

## 2. Reference / annotation

- **Assembly:** hg19 / GRCh37 (`@SQ SN:chr1 … chrX` from STAR alignment).
- **Aligner:** STAR 2.7.5a, two-pass basic, sorted-by-coordinate, dup-marked
  with `sambamba markdup` (per `@PG` lines in BAM headers).
- **Annotation used for sashimi:** GENCODE v19 (hg19, `chr`-prefixed) —
  downloaded to `../annotation/gencode.v19.annotation.gtf.gz`.
  Do **NOT** substitute a GRCh38 (hg38) GTF such as
  `Homo_sapiens.GRCh38.109.gtf.gz`: it refers to a different assembly and
  also uses `X` instead of `chrX`, so features won't line up with the BAM
  coordinates.

---

## 3. Sample sheets for ggsashimi

ggsashimi expects a headerless TSV: `sample_id<TAB>bam_path<TAB>color[<TAB>group]`.

| File | Description |
|------|-------------|
| `samples.tsv`  | 3-column TSV (`id<TAB>bam<TAB>color_factor`) with **absolute** BAM paths. Absolute paths are required inside the Docker container — relative paths fail with `No available bam files.` |
| `palette.tsv`  | Palette file passed via `-P`. One color per line, column 1 only (extra tab-separated columns are treated as comments). Line *i* maps to the *i*-th **unique** value of column 3 of `samples.tsv` in order of first appearance. See section 6 for the full mapping. |

### Color scheme

Simplified palette emphasising treatment effect rather than replicate /
emetine differences: all controls near-black, untreated patient red,
ASO-treated patient blue.

| Group                | Factor (col 3 of `samples.tsv`) | Rendered hex (line of `palette.tsv`) | Name        | Semantic                    |
|----------------------|---------------------------------|--------------------------------------|-------------|-----------------------------|
| `AS_015_En` (R1–R3)  | `#B22222`                       | `#CD5C5C` (line 1)                   | indianred   | patient, untreated, −emetine |
| `AS_015_Ep` (R1–R3)  | `#660000`                       | `#8B0000` (line 2)                   | darkred     | patient, untreated, +emetine |
| `control_1_56_*`     | `#000000`                       | `#121212` (line 3)                   | near-black  | control                      |
| `control_1_29_*`     | `#000000`                       | `#121212` (line 3)                   | near-black  | control                      |
| `AS_015_ASO` (R1–R3) | `#000080`                       | `#1E90FF` (line 4)                   | dodgerblue  | patient + ASO                |

> En/Ep tracks within the control groups share a color — the distinction
> is carried by the track label printed on each track. For the patient
> `AS_015` samples, En and Ep are split (indianred vs. darkred) to
> highlight NMD-sensitive transcripts. If you want En vs Ep visually
> differentiated for controls too, give each its own level in column 3
> of `samples.tsv` and add a matching line to `palette.tsv`.

---

## 4. Region of interest

**`chrX:107,811,500–107,815,000`** (hg19, COL4A5, + strand per GENCODE v19).

Features inside the window (GENCODE v19):

| Feature | hg19 coords | Note |
|---|---|---|
| Exon 5 | 107,811,859–107,811,903 | canonical |
| Exon 6 | 107,811,989–107,812,051 | canonical (COL4A5-001 / -006) |
| **Exon 6 extended** | **107,811,989–107,812,265** | **retained-intron tx `ENST00000470339` — pseudoexon / intron-6 inclusion event** |
| Exon 7 | 107,814,643–107,814,696 | canonical |

---

## 5. Tool (ggsashimi)

```
docker pull guigolab/ggsashimi
```

Image: `guigolab/ggsashimi:latest` (linux/amd64; runs on Apple Silicon via
`--platform linux/amd64` emulation).

---

## 6. Analysis commands
### larger per-track height and font
Output: `sashimi_COL4A5_intron6.pdf` (~14" × 53" canvas)

```bash
docker run --rm --platform linux/amd64 -w "$PWD" -v "$PWD":"$PWD" \
  guigolab/ggsashimi \
  -b "$PWD/Targeted_rnaseq/samples.tsv" \
  -c chrX:107811500-107815000 \
  -M 15 \
  -C 3 \
  -P "$PWD/Targeted_rnaseq/palette.tsv" \
  --alpha 0.9 \
  -g "$PWD/annotation/gencode.v19.annotation.gtf.gz" \
  --height 2.5 \
  --width 10 \
  --ann-height 3 \
  --base-size 20 \
  -o sashimi_COL4A5_intron6
```

Notes:
- Pass the TSV and GTF paths as **absolute** (`"$PWD/..."`). Inside the
  container, relative paths like `Targeted_rnaseq/samples.tsv` can fail with
  `FileNotFoundError` depending on how Docker resolves the working directory.
- The BAM paths **inside `samples.tsv`** must also be absolute — ggsashimi
  resolves them relative to the TSV file's directory, not the container CWD.
- `-C 3` tells ggsashimi that column 3 of `samples.tsv` is the color
  **factor**. The hex strings there are treated as **group labels**, not as
  literal colors. The actual plotted colors come from the palette file
  passed via **`-P`** — without `-P` ggsashimi falls back to its hard-coded
  default `"#ff0000", "#00ff00", "#0000ff", "#000000"` (see `read_palette`
  in `ggsashimi.py`) and your column-3 hexes are ignored.
- `palette.tsv` lists **one color per line** (column 1 only; extra
  tab-separated columns are ignored and useful as comments). Line *i* of
  the palette is matched to the *i*-th unique value of the color-factor
  column in its **order of first appearance in `samples.tsv`** (confirmed
  in `ggsashimi.py` via `levels = list(OrderedDict.fromkeys(d.values()).keys())`
  and `p[levels.index(v)]`). For this project column 3 of `samples.tsv`
  contains 4 unique levels in the order
  `#B22222, #660000, #000000, #000080`, so `palette.tsv` has exactly 4
  lines whose line numbers map to:

  | Line | Factor level | Rendered color | Groups hit |
  |------|--------------|----------------|------------|
  | 1    | `#B22222`    | `#CD5C5C` indianred       | `AS_015_En` (R1–R3) |
  | 2    | `#660000`    | `#8B0000` darkred         | `AS_015_Ep` (R1–R3) |
  | 3    | `#000000`    | `#121212` near-black      | all `control_1_56_*` and `control_1_29_*` |
  | 4    | `#000080`    | `#1E90FF` dodgerblue      | `AS_015_ASO` (R1–R3) |

  If you prefer **identity mapping** (rendered color ≡ the hex you wrote
  in column 3), just copy those 4 hexes verbatim into `palette.tsv` in the
  same order.
- `-o sashimi_COL4A5_intron6` (no extension) → ggsashimi appends `.pdf`,
  producing `sashimi_COL4A5_intron6.pdf`.
- `--platform linux/amd64` is required on Apple Silicon (M-series) Macs
  because the published image is amd64-only.

---

## 7. Parameter glossary

| Flag | Meaning | Used value |
|------|---------|------------|
| `-b` | sample TSV (`id<TAB>bam<TAB>color[<TAB>group]`) | `Targeted_rnaseq/samples.tsv` |
| `-c` | genomic region | `chrX:107811500-107815000` |
| `-g` | GTF annotation (gzipped OK) | `annotation/gencode.v19.annotation.gtf.gz` (hg19) |
| `-M` | min reads to draw a junction arc | `15` |
| `-C` | TSV column (1-based) defining color groups | `3` (hex column — used as a *factor*, not as literal color) |
| `-P` | palette file: one color per line, line *i* → *i*-th unique level of `-C` column in order of first appearance | `Targeted_rnaseq/palette.tsv` |
| `-O` | TSV column (1-based) for **overlay** (aggregating samples per track) | not used in the final runs (would require adding a 4th column with condition labels to `samples.tsv`) |
| `-A` | aggregation function for `-O`: `mean / median / mean_j / median_j` | not used (kept per-replicate tracks) |
| `--alpha` | coverage transparency (0–1) | `0.9` |
| `--shrink` | visually compress long introns | not used |
| `--height` | height per track (inches) | `2.5` (default 2) |
| `--width` | total figure width (inches) | `10` (default 10) |
| `--ann-height` | height of GTF annotation track (inches) | `3` (default 1.5) |
| `--base-size` | base ggplot font size (pt) | `20` (default 14) |
| `-o` | output prefix → writes `<prefix>.pdf` | run-specific |

---

## 8. Caveats

- The htslib warnings `[W::hts_idx_load3] The index file is older than the data file: …bam.bai` are harmless — the BAM files themselves are unchanged, only file mtimes differ. To silence:

  ```bash
  for b in Targeted_rnaseq/bams/*.bam; do samtools index "$b"; done
  ```

- `--platform linux/amd64` is required on Apple Silicon (M-series) Macs because
  the published image is amd64-only. Performance is fine via Rosetta/QEMU
  emulation for these small targeted BAMs.

---

## 9. Outputs

Generated in the repository root:

| File | Notes |
|------|-------|
| `sashimi_COL4A5_intron6.pdf` | Vector output from the command in section 6 
