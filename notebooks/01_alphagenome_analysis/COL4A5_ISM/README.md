# Download the UCSC RefSeq annotation table for hg38
wget "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"
gunzip ncbiRefSeq.txt.gz

# Download hg38 genome sequence FASTA file
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Check if NM_033380.3 is there
grep "NM_033380" ncbiRefSeq.txt

# Run script to generate bed and fasta file containing all intronic sequences 
python3 extract_introns.py

# Running AlphaGenome on all introns from COL4A5 gene 
jupyter execute alphagenome_ISM_COL4A5.ipynb 