## Run RAxML
cp ../core_gene_identification/Output/core_gene_alignment.aln .
trimal -in core_gene_alignment.aln -out core_gene_alignment_trimmed.aln -automated1 -fasta
raxmlHPC-HYBRID-SSE3 -T 25 -s core_gene_alignment_trimmed.aln -N 100 -n coreGeneAlignment -f a -p 12345 -x 12345 -m GTRGAMMA
