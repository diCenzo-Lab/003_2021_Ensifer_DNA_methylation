## Run RAxML
cp ../core_gene_identification/Output/core_gene_alignment.aln .
raxmlHPC-HYBRID-SSE3 -T 25 -s core_gene_alignment.aln -N autoMRE -n coreGeneAlignment -f a -p 12345 -x 12345 -m GTRGAMMA
