#!/usr/bin/env python3
"""
Extract intronic sequences for any RefSeq transcript from a genome assembly.

Given an NCBI RefSeq transcript ID (e.g. NM_033380.3), this script parses the
UCSC ncbiRefSeq table, derives the intron coordinates, and writes per-intron
BED/TSV/FASTA outputs. Output file names embed both the gene name and the
transcript ID, e.g. ``COL4A5_NM_033380.3_introns.bed``.

Inputs:
    - ncbiRefSeq.txt   (UCSC RefSeq table)
    - genome FASTA     (e.g. hg38.fa) for sequence extraction (optional)

Outputs (in --outdir):
    - <GENE>_<TRANSCRIPT>_introns.bed   (coordinates + length)
    - <GENE>_<TRANSCRIPT>_introns.tsv   (detailed table)
    - <GENE>_<TRANSCRIPT>_introns.fa    (sequences, if genome FASTA available)

Examples:
    python extract_introns.py NM_033380.3
    python extract_introns.py NM_000546 --refseq ncbiRefSeq.txt --genome hg38.fa
    python extract_introns.py NM_004006.3 --outdir results/ --gene DMD
"""

import argparse
import os
import re
import subprocess
import sys


# ─────────────────────────────────────────────
# Parse RefSeq table
# ─────────────────────────────────────────────
def parse_refseq_table(filename, transcript_id):
    """
    Parse UCSC ncbiRefSeq table and return matching transcript records.

    Column layout (0-based index):
    0:  bin
    1:  name        ← transcript ID (NM_...)
    2:  chrom
    3:  strand
    4:  txStart     ← 0-based transcript start
    5:  txEnd       ← transcript end
    6:  cdsStart
    7:  cdsEnd
    8:  exonCount
    9:  exonStarts  ← comma-separated, 0-based
    10: exonEnds    ← comma-separated
    12: name2       ← gene name
    """
    found = []
    base_id = transcript_id.split(".")[0]

    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            name = fields[1]

            # Match exact (with version) or base accession (without version)
            if name == transcript_id or name.split(".")[0] == base_id:
                exon_starts = [int(x) for x in fields[9].rstrip(",").split(",")]
                exon_ends   = [int(x) for x in fields[10].rstrip(",").split(",")]
                gene_name   = fields[12] if len(fields) > 12 else "unknown"

                found.append({
                    "name"        : name,
                    "gene_name"   : gene_name,
                    "chrom"       : fields[2],
                    "strand"      : fields[3],
                    "tx_start"    : int(fields[4]),   # 0-based
                    "tx_end"      : int(fields[5]),
                    "cds_start"   : int(fields[6]),
                    "cds_end"     : int(fields[7]),
                    "exon_count"  : int(fields[8]),
                    "exon_starts" : exon_starts,
                    "exon_ends"   : exon_ends,
                })

    return found


def select_transcript(transcripts, transcript_id):
    """
    Choose a single transcript record from possibly several matches.

    Preference order:
    1. Exact match on the full transcript ID (including version).
    2. Otherwise the first record (base-accession match).

    Alternate/patch contigs (chrom containing '_') are de-prioritised so the
    primary assembly placement is preferred when both exist.
    """
    if not transcripts:
        return None

    def is_primary(t):
        return "_" not in t["chrom"]

    exact = [t for t in transcripts if t["name"] == transcript_id]
    pool = exact if exact else transcripts

    primary = [t for t in pool if is_primary(t)]
    chosen_pool = primary if primary else pool

    if len(chosen_pool) > 1:
        locations = ", ".join(
            f"{t['name']}@{t['chrom']}:{t['tx_start']}-{t['tx_end']}" for t in chosen_pool
        )
        print(f"WARNING: multiple matches found, using the first: {locations}")

    return chosen_pool[0]


# ─────────────────────────────────────────────
# Derive intron coordinates
# ─────────────────────────────────────────────
def get_introns(transcript):
    """
    Introns are gaps BETWEEN consecutive exons.

    Exon structure (0-based BED coords):
    Exon1: [exon_starts[0], exon_ends[0]]
    Intron1: [exon_ends[0], exon_starts[1]]
    Exon2: [exon_starts[1], exon_ends[1]]
    ...
    """
    exon_starts = transcript["exon_starts"]
    exon_ends   = transcript["exon_ends"]
    strand      = transcript["strand"]
    chrom       = transcript["chrom"]
    introns     = []

    n_introns = len(exon_ends) - 1  # Number of introns = exons - 1

    for i in range(n_introns):
        intron_start = exon_ends[i]        # 0-based start
        intron_end   = exon_starts[i + 1]  # end
        intron_len   = intron_end - intron_start

        # Intron number relative to strand
        if strand == "+":
            intron_num = i + 1
        else:
            # For minus strand, reverse numbering
            intron_num = n_introns - i

        introns.append({
            "chrom"       : chrom,
            "start"       : intron_start,   # 0-based (BED)
            "end"         : intron_end,
            "name"        : f"intron_{intron_num}",
            "length"      : intron_len,
            "strand"      : strand,
            "intron_index": i + 1           # raw index for reference
        })

    return introns


# ─────────────────────────────────────────────
# Output naming
# ─────────────────────────────────────────────
def sanitize(text):
    """Make a string safe for use in a file name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")


def build_output_paths(outdir, gene_name, transcript_name, prefix=None):
    """Construct the BED/TSV/FASTA output paths embedding gene + transcript."""
    if prefix is None:
        prefix = f"{sanitize(gene_name)}_{sanitize(transcript_name)}"
    base = os.path.join(outdir, f"{prefix}_introns")
    return {
        "prefix": prefix,
        "bed": f"{base}.bed",
        "tsv": f"{base}.tsv",
        "fasta": f"{base}.fa",
    }


# ─────────────────────────────────────────────
# Write output files
# ─────────────────────────────────────────────
def write_bed(introns, filename):
    """Write introns to BED6 format"""
    with open(filename, "w") as out:
        for intron in introns:
            out.write(
                f"{intron['chrom']}\t"
                f"{intron['start']}\t"
                f"{intron['end']}\t"
                f"{intron['name']}\t"
                f"{intron['length']}\t"
                f"{intron['strand']}\n"
            )
    print(f"BED file written: {filename}")


def write_tsv(introns, transcript, filename):
    """Write detailed intron table"""
    with open(filename, "w") as out:
        # Header
        out.write(
            "intron_name\tchrom\tstart_hg38\tend_hg38\t"
            "length_bp\tstrand\ttranscript\tgene\n"
        )
        for intron in introns:
            out.write(
                f"{intron['name']}\t"
                f"{intron['chrom']}\t"
                f"{intron['start']}\t"      # 0-based BED
                f"{intron['end']}\t"
                f"{intron['length']}\t"
                f"{intron['strand']}\t"
                f"{transcript['name']}\t"
                f"{transcript['gene_name']}\n"
            )
    print(f"TSV file written: {filename}")


def extract_sequences(bed_file, genome_fasta, output_fasta):
    """Use bedtools getfasta to extract intronic sequences"""

    # Check bedtools is available
    result = subprocess.run(["which", "bedtools"], capture_output=True, text=True)
    if result.returncode != 0:
        print("bedtools not found — skipping FASTA extraction")
        return

    # Index FASTA if needed
    if not os.path.exists(f"{genome_fasta}.fai"):
        print("Indexing FASTA...")
        subprocess.run(["samtools", "faidx", genome_fasta])

    # Extract sequences
    cmd = [
        "bedtools", "getfasta",
        "-fi", genome_fasta,
        "-bed", bed_file,
        "-s",              # strand-aware (reverse complement for minus strand)
        "-name",           # use name column as FASTA header
        "-fo", output_fasta
    ]

    subprocess.run(cmd, check=True)
    print(f"FASTA written: {output_fasta}")


# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────
def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract intron coordinates/sequences for a RefSeq transcript.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "transcript",
        help="NCBI RefSeq transcript ID, with or without version (e.g. NM_033380.3).",
    )
    parser.add_argument(
        "--refseq", default="ncbiRefSeq.txt",
        help="Path to the UCSC ncbiRefSeq table.",
    )
    parser.add_argument(
        "--genome", default="hg38.fa",
        help="Path to the genome FASTA used for sequence extraction.",
    )
    parser.add_argument(
        "--outdir", default=".",
        help="Directory where output files are written.",
    )
    parser.add_argument(
        "--gene", default=None,
        help="Override the gene name used in output file names "
             "(default: gene name from the RefSeq table).",
    )
    parser.add_argument(
        "--prefix", default=None,
        help="Override the full output file prefix "
             "(default: <GENE>_<TRANSCRIPT>).",
    )
    return parser.parse_args(argv)


# ─────────────────────────────────────────────
# Main function
# ─────────────────────────────────────────────
def main(argv=None):
    args = parse_args(argv)
    transcript_id = args.transcript

    print(f"\n{'='*50}")
    print(f"Extracting introns for: {transcript_id}")
    print(f"{'='*50}\n")

    if not os.path.exists(args.refseq):
        print(f"ERROR: RefSeq table not found: {args.refseq}")
        return 1

    # Parse RefSeq table
    transcripts = parse_refseq_table(args.refseq, transcript_id)

    if not transcripts:
        print(f"ERROR: {transcript_id} not found in {args.refseq}")
        print("Tip: try without the version number (e.g. NM_033380)")
        return 1

    transcript = select_transcript(transcripts, transcript_id)

    # Allow the user to override the gene name used for file naming
    gene_name = args.gene if args.gene else transcript["gene_name"]

    os.makedirs(args.outdir, exist_ok=True)
    paths = build_output_paths(args.outdir, gene_name, transcript["name"], args.prefix)

    print(f"Transcript   : {transcript['name']}")
    print(f"Gene         : {transcript['gene_name']}")
    print(f"Location     : {transcript['chrom']}:"
          f"{transcript['tx_start']}-{transcript['tx_end']}")
    print(f"Strand       : {transcript['strand']}")
    print(f"Exon count   : {transcript['exon_count']}")
    print(f"Intron count : {transcript['exon_count'] - 1}")
    print(f"Output prefix: {paths['prefix']}\n")

    # Get introns
    introns = get_introns(transcript)

    if not introns:
        print("No introns found (single-exon transcript). Nothing to write.")
        return 0

    # Print summary table to screen
    print(f"{'Intron':<12} {'Chrom':<6} {'Start (hg38)':<15} "
          f"{'End (hg38)':<15} {'Length (bp)':<12} {'Strand'}")
    print("-" * 65)

    for intron in introns:
        print(
            f"{intron['name']:<12} "
            f"{intron['chrom']:<6} "
            f"{intron['start']:<15} "
            f"{intron['end']:<15} "
            f"{intron['length']:<12} "
            f"{intron['strand']}"
        )

    total_intronic = sum(i['length'] for i in introns)
    print(f"\nTotal intronic sequence: {total_intronic:,} bp")
    print(f"Number of introns:       {len(introns)}")

    # Write files
    print()
    write_bed(introns, paths["bed"])
    write_tsv(introns, transcript, paths["tsv"])

    # Extract sequences if genome FASTA exists
    if os.path.exists(args.genome):
        extract_sequences(paths["bed"], args.genome, paths["fasta"])
    else:
        print(f"\nWARNING: {args.genome} not found")
        print("Skipping sequence extraction")
        print("Run manually:\n  bedtools getfasta "
              f"-fi {args.genome} -bed {paths['bed']} -s -name -fo {paths['fasta']}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
