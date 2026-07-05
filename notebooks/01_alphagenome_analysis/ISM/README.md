# Full-gene intron ISM (COL4A5, and any RefSeq transcript)

This folder builds a per-intron saturation-mutagenesis (ISM) workflow for **all
introns of a gene** on hg38. Although the manuscript uses `COL4A5`
(`NM_033380.3`), both the extraction script and the scoring notebook are
gene-agnostic — set the transcript ID / gene name and rerun.

## Files

- `extract_introns.py`
  - Command-line tool that parses the UCSC `ncbiRefSeq` table for a given
    transcript, derives intron coordinates, and writes per-intron
    `BED`/`TSV`/`FASTA` files. Output names embed gene + transcript, e.g.
    `COL4A5_NM_033380.3_introns.{bed,tsv,fa}`.
- `alphagenome_ISM_batch.ipynb`
  - Gene-agnostic **batch** ISM notebook. For each intron it builds a VCF of
    every single-nucleotide substitution (`length_bp × 3` rows), scores them in
    parallel with the AlphaGenome splice scorers (`SPLICE_JUNCTIONS`,
    `SPLICE_SITE_USAGE`, `SPLICE_SITES`), filters to the gene + tissue ontology,
    and writes one tidy CSV per intron plus a combined CSV.
  - Set `GENE_NAME` / `TRANSCRIPT_ID` / `ONTOLOGY_CURIE` in the **Configuration**
    cell (defaults: `COL4A5`, `NM_033380.3`, `UBERON:0002113` = kidney).
  - Outputs go to `../alphagenome_ISM_<GENE>/` (gitignored), with one
    `intron_*/` subfolder each containing `<intron>_variants.vcf` and
    `csv/alphagenome_scores_<intron>.csv`. This is the faster path based on the
    [AlphaGenome batch variant scoring](https://www.alphagenomedocs.com/colabs/batch_variant_scoring.html)
    tutorial; introns whose output folder already exists are skipped unless
    `overwrite=True`.

## Prepare inputs

```bash
# 1. Download the UCSC RefSeq annotation table for hg38
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
gunzip ncbiRefSeq.txt.gz

# 2. Download the hg38 genome FASTA (needed for intron sequence extraction)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# 3. (optional) confirm the transcript is present
grep "NM_033380" ncbiRefSeq.txt

# 4. Generate BED / TSV / FASTA of all intronic sequences for the transcript
python3 extract_introns.py NM_033380.3
#   → COL4A5_NM_033380.3_introns.{bed,tsv,fa}
#
#   Other examples:
#     python3 extract_introns.py NM_000091.5 --gene COL4A3
#     python3 extract_introns.py NM_004006.3 --outdir results/ --gene DMD
#   Run `python3 extract_introns.py --help` for all options
#   (--refseq, --genome, --outdir, --gene, --prefix).
```

`extract_introns.py` uses `bedtools getfasta` (strand-aware) for sequence
extraction; if `bedtools`/`samtools` or the genome FASTA are missing it still
writes the BED/TSV and prints the manual command to run later.

## Run the ISM scoring

Open `alphagenome_ISM_batch.ipynb`, set the configuration cell (and your
AlphaGenome API key), then run all cells:

```bash
jupyter execute alphagenome_ISM_batch.ipynb
```

Downstream summaries and cross-intron hotspot comparisons are produced by
`../../02_postprocessing_plots_summary/alphagenome_analysis_notebook_ISM_all_introns.ipynb`.
