#!usr/bin/perl
use 5.010;

# Get files
system("mkdir Sinorhizobium"); # Make directory to hold the genome fna files
system("mkdir Sinorhizobium2"); # Make directory to hold the CDS fna files
system("perl downloadGenomes.pl genomeList_sino.txt"); # Get the genome fna files
system("perl downloadGenomes2.pl genomeList_sino.txt"); # Get the CDS fna files

# make an array of the species names
$species_list = 'genomeList_sino.txt';
open($species,'<',$species_list);
while(<$species>) {
	@line = split("\t",$_); # split the input line into an array based on commas
	push(@genomes,@line[0]); # make an array of the species names
}
close($species);

# Get GANTC count per genome, and also record genome size and GC content
foreach $i (@genomes) {
	$genome2 = 'Sinorhizobium/' . $i . '.fna';
	system("printf $i >> output_wholeGenome.txt");
	system("perl countGANTC.pl $genome2 >> output_wholeGenome.txt");
}

# Get GANTC count for coding regions per genome, and also record genome size and GC content
foreach $i (@genomes) {
	$genome2 = 'Sinorhizobium2/' . $i . '.fna';
	system("printf $i >> output_CDS.txt");
	system("perl countGANTC.pl $genome2 >> output_CDS.txt");
}
