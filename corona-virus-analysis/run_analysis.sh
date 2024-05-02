#!/usr/bin/env bash

# store SRR ID and other variables
SRR="SRR10971381"

# retrieve SRA metadata
# bio search $SRR > $SRR.metadata

# create directory for reads
mkdir -p reads

# specify a subset size
# N=$1

# obtain SRA reads
echo "Downloading $N reads..."
fastq-dump $SRR --split-files --origfmt --outdir reads -v

# display summary statistics for reads
seqkit stats reads/$SRR*

# remove adaptor sequences
F1="reads/${SRR}_1.fastq"
F2="reads/${SRR}_2.fastq"

echo "Trimming adaptor sequences..."
trimmomatic PE $F1 $F2 -baseout reads/read.fq SLIDINGWINDOW:4:30

# generate quality control reports
# fastqc reads/*.fq

# assemble reads using megahit
R1="reads/read_1P.fq"
R2="reads/read_2P.fq"

rm -rf out	

echo "Assembling contigs..."
megahit -1 $R1 -2 $R2 -o out

# generate assembly statistics
echo "ASSEMBLY STATISTICS:"
seqkit stats out/final.contigs.fa

# create BLAST database from generated contigs
makeblastdb -out out/contigs -in out/final.contigs.fa -dbtype nucl -parse_seqids

# query the database to retrieve largest contigs
blastdbcmd -db out/contigs -entry all -outfmt "%l %a" | sort -rn | head

# extract the sequence that corresponds to the largest contig
CONTIG_ID=$(blastdbcmd -db out/contigs -entry all -outfmt "%l %a" | sort -rn | head -n 1 | awk '{ print $2 }')
blastdbcmd -db out/contigs -entry $CONTIG_ID > genome.fa

# run a BLAST query on the longset contig
blastn -db db/ref_viruses_rep_genomes -query genome.fa -outfmt "6 pident length sacc stitle" > blast.txt
cat blast.txt

