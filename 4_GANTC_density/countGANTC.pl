#!usr/bin/perl
use 5.010;

# Set variables
$count = 0; # A variable meant to allow the first line to be treated differently from the rest
$GANTC = 0; # The number of GANTC sites
$GC = 0; # The number of G or C nucleotides
$nt = 0; # The total number of nucleotides
$temp_file = 'temporary.txt'; # A temporary file

# Store the sequences in a temp file, with each replicon/CDS on a separate line and with the full replicon/CDS on one line
open($out, '>', $temp_file);
while(<>) {
  chomp;
  if(/>/) {
    $count++;
    print $out ("\n") unless($count == 1); # If a new replicon/CDS, move to a new line in the temp file
  }
  else {
    print $out ("$_"); # Print the sequence to the temp file
  }
}
close($out);

# Determine GANTC count, GC content, and genome length (or total CDS length)
open($in, '<', $temp_file);
while(<$in>) {
  @line = split('', $_); # Store each nt as a separate entry in an array
  $length = scalar(@line); # Determine length of the array
  for($n=0; $n < $length - 4; $n++) { # Loop through the array (i.e., the sequence), and count all occurances of the GANTC motif
    $test = @line[$n] . @line[$n+1] . @line[$n+2] . @line[$n+3] . @line[$n+4];
    if($test eq 'GAATC' || $test eq 'GAGTC' || $test eq 'GATTC' || $test eq 'GACTC') {
      $GANTC++;
    }
  }
  for($n=0; $n < $length; $n++) { # Loop through the array (i.e., the sequence), and count all G or C nts and the total number of nts
    $nt++;
    if(@line[$n] eq 'G' || @line[$n] eq 'C') {
      $GC++;
    }
  }
}
close($in);

# Record values
$density = 1000 * $GANTC / $nt; # Determine the frequency of GANTC sites, as number of GANTC sites per kb
$GC_content = 100 * $GC / $nt; # Determine the GC content
say("\t$GANTC\t$nt\t$density\t$GC_content"); # Record the total number of GANTC sites, total genome/CDS length, the frequency of GANTC sites, and the GC content
unlink($temp_file);
