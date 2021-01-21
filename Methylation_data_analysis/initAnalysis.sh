#!usr/bin/bash

# Help message
if [ "$1" == "-h" ]; then
  echo "
  Usage: `basename $0` motifs.csv motifs.gff genome.gff genome.fna clean

  Where:
    motifs.csv: A table containing summary information for each of the methylated motifs detected in the genome.
    motifs.gff: A gff file providing information on each of the methylated sites in the genome.
    genome.gff: The genome GFF file.
    genome.fna: The genome sequence fasta file.
    clean: A value of 1 or 0 to indicate if intermediate files should be remove (1) or kept (0).
    "
  exit 0
fi
if [ "$1" == "-help" ]; then
  echo "
  Usage: `basename $0` motifs.csv motifs.gff genome.gff genome.fna clean

  Where:
    motifs.csv: A table containing summary information for each of the methylated motifs detected in the genome.
    motifs.gff: A gff file providing information on each of the methylated sites in the genome.
    genome.gff: The genome GFF file.
    genome.fna: The genome sequence fasta file.
    clean: A value of 1 or 0 to indicate if intermediate files should be remove (1) or kept (0).
    "
    exit 0
fi

# Add folders if they don't exist
if [ ! -d "Intermediate" ]; then
  mkdir Intermediate
fi
if [ ! -d "Output" ]; then
  mkdir Output
fi

# Separate methylation data for each replicon
perl Scripts/extractReplicons.pl $1 $2 $3

# Calculate methylation extent across genome
perl Scripts/extractMethylationSkew.pl $1 $2 $3

# Calculate GC skew across genome
perl Scripts/extractGcSkew.pl $1 $2 $4

# Calculate methylation extent of each gene
perl Scripts/extractGeneMethylation.pl $1 $2 $3

# Calculate methylation extent of the upstream region of each gene
perl Scripts/extractPromoterMethylation.pl $1 $2 $3 100
perl Scripts/extractPromoterMethylation.pl $1 $2 $3 125
perl Scripts/extractPromoterMethylation.pl $1 $2 $3 250

# Clean
if [ "$5" == 1 ]; then
  rm -r Intermediate
fi
