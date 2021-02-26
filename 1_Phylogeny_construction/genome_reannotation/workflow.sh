## Prepare necessary directories
mkdir Genomes/
mkdir Output/
mkdir intermediaryFiles/

## Get the genomes list
mv refseqGenomeInformation.txt intermediaryFiles/refseqGenomeInformation.txt # Move the file
perl Scripts/parseGenomeList.pl intermediaryFiles/refseqGenomeInformation.txt # Parse the NCBI Genome information
mv genomeList.txt intermediaryFiles/ # Move the file
perl Scripts/downloadGenomes.pl intermediaryFiles/genomeList.txt # Get the genomes

## Reannotate the genomes
perl Scripts/runProkka.pl intermediaryFiles/genomeList.txt # Run prokka to annotate the genomes
