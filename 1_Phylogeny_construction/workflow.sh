## Reannotate the genomes
cd genome_reannotation # Move into the appropriate directory
sh workflow.sh # Run the analysis
cd .. # Move back into the main directory

## Identify the core genes
cd core_gene_identification # Move into the appropriate directory
sh workflow.sh # Run the analysis
cd .. # Move back into the main directory

## Make the phylogeny
cd core_gene_phylogeny # Move into the appropriate directory
sh workflow.sh # Run the analysis
cd .. # Move back into the main directory

## Identify symbiotic genes
cd symbiotic_genes # Move into the appropriate directory
sh workflow.sh # Run the analysis
cd .. # Move back into the main directory

