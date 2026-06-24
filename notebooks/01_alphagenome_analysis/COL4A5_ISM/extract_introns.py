#!/usr/bin/env python3
"""
Extract intronic sequences for NM_033380.3 (COL4A5) from hg38
Input:  - ncbiRefSeq.txt (UCSC RefSeq table)
        - hg38.fa (genome FASTA)
Output: - NM_033380_introns.bed      (coordinates + length)
        - NM_033380_introns.tsv      (detailed table)
        - NM_033380_introns.fa       (sequences)
"""

import os
import re
import subprocess

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
TRANSCRIPT_ID  = "NM_033380.3"
REFSEQ_TABLE   = "ncbiRefSeq.txt"
GENOME_FASTA   = "hg38.fa"
OUTPUT_BED     = "NM_033380_introns.bed"
OUTPUT_TSV     = "NM_033380_introns.tsv"
OUTPUT_FASTA   = "NM_033380_introns.fa"

# ─────────────────────────────────────────────
# Parse RefSeq table
# ─────────────────────────────────────────────
def parse_refseq_table(filename, transcript_id):
    """
    Parse UCSC ncbiRefSeq table and return transcript info.
    
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
    
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue
            
            name = fields[1]
            
            # Match exact or base accession (without version)
            base_id = transcript_id.split(".")[0]
            if name == transcript_id or name == base_id:
                
                chrom       = fields[2]
                strand      = fields[3]
                tx_start    = int(fields[4])   # 0-based
                tx_end      = int(fields[5])
                cds_start   = int(fields[6])
                cds_end     = int(fields[7])
                exon_count  = int(fields[8])
                
                # Parse exon starts and ends (comma-separated, trailing comma)
                exon_starts = [int(x) for x in fields[9].rstrip(",").split(",")]
                exon_ends   = [int(x) for x in fields[10].rstrip(",").split(",")]
                gene_name   = fields[12] if len(fields) > 12 else "unknown"
                
                found.append({
                    "name"        : name,
                    "gene_name"   : gene_name,
                    "chrom"       : chrom,
                    "strand"      : strand,
                    "tx_start"    : tx_start,
                    "tx_end"      : tx_end,
                    "cds_start"   : cds_start,
                    "cds_end"     : cds_end,
                    "exon_count"  : exon_count,
                    "exon_starts" : exon_starts,
                    "exon_ends"   : exon_ends,
                })
    
    return found

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
# Main function
# ─────────────────────────────────────────────
def main():
    
    print(f"\n{'='*50}")
    print(f"Extracting introns for: {TRANSCRIPT_ID}")
    print(f"{'='*50}\n")
    
    # Parse RefSeq table
    transcripts = parse_refseq_table(REFSEQ_TABLE, TRANSCRIPT_ID)
    
    if not transcripts:
        print(f"ERROR: {TRANSCRIPT_ID} not found in {REFSEQ_TABLE}")
        print("Tip: try without version number (NM_033380)")
        return
    
    # Use first match (should only be one)
    transcript = transcripts[0]
    
    print(f"Transcript   : {transcript['name']}")
    print(f"Gene         : {transcript['gene_name']}")
    print(f"Location     : {transcript['chrom']}:"
          f"{transcript['tx_start']}-{transcript['tx_end']}")
    print(f"Strand       : {transcript['strand']}")
    print(f"Exon count   : {transcript['exon_count']}")
    print(f"Intron count : {transcript['exon_count'] - 1}\n")
    
    # Get introns
    introns = get_introns(transcript)
    
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
    write_bed(introns, OUTPUT_BED)
    write_tsv(introns, transcript, OUTPUT_TSV)
    
    # Extract sequences if genome FASTA exists
    if os.path.exists(GENOME_FASTA):
        extract_sequences(OUTPUT_BED, GENOME_FASTA, OUTPUT_FASTA)
    else:
        print(f"\nWARNING: {GENOME_FASTA} not found")
        print("Skipping sequence extraction")
        print(f"Run manually:\n  bedtools getfasta -fi hg38.fa -bed {OUTPUT_BED} -s -name -fo {OUTPUT_FASTA}")

if __name__ == "__main__":
    main()
