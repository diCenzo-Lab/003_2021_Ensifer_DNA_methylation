cat *.faa > combinedProteomes.fasta
cd-hit -i combinedProteome.fasta -o output -c 0.7 -G 0 -M 50000 -T 20 -n 4 -d 0 -aL 0.8
