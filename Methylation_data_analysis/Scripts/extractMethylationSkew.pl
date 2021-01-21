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

# Extract replicons
open($in, '<', $gffFiles);
while(<$in>) {
  chomp;
  if(/sequence-region/) {
    @line = split(' ', $_);
    push(@replicons, @line[1]);
  }
}

# Bin into 10kb bins with a 10 kb sliding window
foreach $i (@motifs) {
  $input = $gffFilesOutput . '_' . $i . '.txt';
  $gffFilesOutput2 = $gffFilesOutput;
  $gffFilesOutput2 =~ s/Intermediate/Output/;
  $output_temp = $gffFilesOutput2 . '_' . $i;
  foreach $j (@replicons) {
    @data = ();
    $output = $output_temp . '_' . $j . '_10kb.txt';
    open($out, '>', $output);
    open($in, '<', $input);
    while(<$in>) {
      if(/\#\#/) {
        if(/$j/) {
          @line = split(' ', $_);
          $size = @line[3];
        }
      }
      else {
          if(/$j/) {
          @line = split("\t", $_);
          @line2 = split(';', @line[8]);
          $methExtent = @line2[7];
          $methExtent =~ s/frac\=//;
          $position = @line[3];
          @data[$position] = $methExtent;
        }
      }
    }
    for($n = 1; $n <= $size; $n = $n + 10000) {
      $methSum = 0;
      $count = 0;
      $min2 = $n - 9999;
      if($min2 < 1) {
        $min = $size + $min2;
      }
      else {
        $min = $min2;
      }
      if($min2 < 1) {
        for($i = 1; $i <= $n; $i++) {
          $methSum = $methSum + @data[$i];
          $test = @data[$i];
          if($test ne "") {
            $count++;
          }
        }
        for($i = $min; $i <= $size; $i++) {
          $methSum = $methSum + @data[$i];
          $test = @data[$i];
          if($test ne "") {
            $count++;
          }
        }
      }
      if($min2 >= 1) {
        for($i = $min; $i <= $n; $i++) {
          $methSum = $methSum + @data[$i];
          $test = @data[$i];
          if($test ne "") {
            $count++;
          }
        }
      }
      if($count == 0) {
        say $out ("$min\t$n\t\t");
      }
      else {
        $avgMeth = $methSum / $count;
        say $out ("$min\t$n\t$avgMeth\t$count");
      }
    }
    $methSum = 0;
    $count = 0;
    for($i = $size - 9999; $i <= $size; $i++) {
      $methSum = $methSum + @data[$i];
      $test = @data[$i];
      if($test ne "") {
        $count++;
      }
    }
    if($count == 0) {
      $temp = $size - 9999;
      say $out ("$temp\t$size\t\t");
    }
    else {
      $avgMeth = $methSum / $count;
      $min = $size - 9999;
      say $out ("$min\t$size\t$avgMeth\t$count");
    }
    close($in);
    close($out);
    $content = '';
    open($in, '<', $output);
    while(<$in>) {
      $content = $content . $_;
    }
    close($in);
    open($out, '>', $output);
    say $out ("Start_nt\tEnd_nt\tExtent_methylation\tMotif_count");
    print $out ($content);
    close($out);
  }
}
