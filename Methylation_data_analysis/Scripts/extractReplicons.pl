#!usr/bin/perl
use 5.010;

# Get file names
$motifFile = @ARGV[0];
$gffFiles = @ARGV[1];
$genomeFile = @ARGV[2];

# Modify gffFiles name in preparation for output files
$gffFilesOutput = $gffFiles;
$gffFilesOutput =~ s/_modified.gff//;
$gffFilesOutput =~ s/Data/Intermediate/;

# Extract motif types
$n = 0;
open($in, '<', $motifFile);
while(<$in>) {
  $_ =~ s/\"//g;
  @line = split(',', $_);
  if(@line[0] eq 'motifString') {
  }
  elsif(@line[0] eq 'C') {
    @motifs[$n] = 'modified_base';
    @motifsType[$n] = 'modified_base';
    $n++;
  }
  else {
    @motifs[$n] = @line[0];
    @motifsType[$n] = @line[2];
    $n++;
  }
}
close($in);

# Extract the relevant data for each motif
$n = 0;
foreach $i (@motifs) {
  $testMotif = 'motif=' . $i;
  $output = $gffFilesOutput . '_' . $i . '.txt';
  open($out, '>', $output);
  open($in, '<', $gffFiles);
  while(<$in>) {
    if(/sequence-region/) {
      print $out ($_);
    }
    if($i eq 'modified_base') {
      if(/modified_base/) {
        print $out ($_);
      }
    }
    else {
      @line = split("\t", $_);
      if(@line[2] eq @motifsType[$n]) {
        if(/$testMotif/) {
          print $out ($_);
        }
      }
    }
  }
  close($out);
  close($in);
  $n++;
}
