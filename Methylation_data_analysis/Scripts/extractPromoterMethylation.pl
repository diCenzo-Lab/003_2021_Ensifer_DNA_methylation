#!usr/bin/perl
use 5.010;

# Get file names
$motifFile = @ARGV[0];
$gffFiles = @ARGV[1];
$genomeFile = @ARGV[2];
$promoter = @ARGV[3];

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

# Extract promoter regions
undef(@geneNames);
undef(@geneStart);
undef(@geneEnd);
undef(@geneReplicon);
open($in, '<', $genomeFile);
while(<$in>) {
  chomp;
  @line = split("\t", $_);
  if(@line[2] eq 'gene') {
    @line2 = split('locus_tag=', @line[8]);
    @line3 = split(';', @line2[1]);
    $name = @line3[0];
  }
  if(@line[2] eq 'CDS') {
    push(@geneNames, $name);
    push(@geneReplicon, @line[0]);
    if(@line[6] eq '+') {
      push(@geneStart, @line[3] - $promoter);
      push(@geneEnd, @line[3] - 1);
    }
    else {
      push(@geneStart, @line[4] + 1);
      push(@geneEnd, @line[4] + $promoter);
    }
  }
}
close($in);

# Extract gene coding regions
open($in, '<', $genomeFile);
while(<$in>) {
  chomp;
  @line = split("\t", $_);
  if(@line[2] eq 'gene') {
    @line2 = split('locus_tag=', @line[8]);
    @line3 = split(';', @line2[1]);
    $name = @line3[0];
  }
  if(@line[2] eq 'CDS') {
    push(@geneNames, $name);
    push(@geneStart, @line[3]);
    push(@geneEnd, @line[4]);
    push(@geneReplicon, @line[0]);

  }
}
close($in);

# Bin by promoter
foreach $i (@motifs) {
  $n = 0;
  $lineNumber = 0;
  $k = 0;
  $input = $gffFilesOutput . '_' . $i . '.txt';
  $gffFilesOutput2 = $gffFilesOutput;
  $gffFilesOutput2 =~ s/Intermediate/Output/;
  $output_temp = $gffFilesOutput2 . '_' . $i;
  foreach $j (@replicons) {
    $testVariable2 = 0;
    $test = 0;
    $testVariable = 0;
    $count = 0; # Number of methylated sites
    $avgMeth = 0; # Average extent of methylation
    $output = $output_temp . '_' . $j . '_promoter' . '_' . $promoter . '.txt';
    @temp = split('_', $j);
    @temp2 = split('\.', @geneReplicon[$n]);
    if(@temp[0] ne @temp2[0]) {
      do {
        $k++;
        $n = $k;
        @temp = split('_', $j);
        @temp2 = split('\.', @geneReplicon[$n]);
      }until(@temp[0] eq @temp2[0])
    }
    open($out, '>', $output);
    open($in, '<', $input);
    while(<$in>) {
      if(/\#\#/) {
      }
      else {
      if($. >= $lineNumber) {
          if(/$j/) {
            @line = split("\t", $_);
            @line2 = split(';', @line[8]);
            $methExtent = @line2[7];
            $methExtent =~ s/frac\=//;
            for($n = $k; $n <= scalar @geneStart - 1; $n++) {
              if(@line[3] >= @geneStart[$n] && @line[3] <= @geneEnd[$n]) {
                $count++;
                $avgMeth = $avgMeth + $methExtent;
                last;
              }
              elsif(@line[3] < @geneStart[$n]) {
                last;
              }
              elsif(@line[3] > @geneEnd[$n]) {
                if($count > 0) {
                  $averageMeth = 100 * $avgMeth / $count;
                  $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
                  $sitesPerGene = 1000 * $count / $geneLength;
                }
                else {
                  $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
                }
                if(/@geneReplicon[$n]/) {
                  $testVariable2 = 0;
                }
                else {
                  $testVariable2 = 1;
                }
                if($testVariable2 == 0) {
                  if($count > 0) {
                    say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t$count\t$sitesPerGene\t$averageMeth");
                  }
                  else {
                    say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t0\t0\tNA");
                  }
                  $k++;
                  $count = 0;
                  $avgMeth = 0;
                }
              }
            }
          }
          else {
            if($test == 0) {
              if($count > 0) {
                $averageMeth = 100 * $avgMeth / $count;
                $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
                $sitesPerGene = 1000 * $count / $geneLength;
                if($testVariable2 == 0) {
                  say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t$count\t$sitesPerGene\t$averageMeth");
                }
              }
              elsif($count == 0) {
                $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
                if($j eq @geneReplicon[$n]) {
                  if($testVariable2 == 0) {
                    say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t0\t0\tNA");
                  }
                }
              }
              $k++;
              $test = 1;
              $lineNumber = $.;
              last;
            }
          }
        }
      }
    }
    @temp = split('_', $j);
    @temp2 = split('\.', @geneReplicon[$n]);
    if(@temp[0] ne @temp2[0]) {
      $k = $k - 1;
    }
    if($count > 0) {
      $averageMeth = 100 * $avgMeth / $count;
      $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
      $sitesPerGene = 1000 * $count / $geneLength;
      if($testVariable2 == 0) {
        say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t$count\t$sitesPerGene\t$averageMeth");
      }
    }
    elsif($count == 0) {
      $geneLength = @geneEnd[$n] - @geneStart[$n] + 1;
      if($j eq @geneReplicon[$n]) {
        if($testVariable2 == 0) {
          say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t0\t0\tNA");
        }
      }
    }
    @temp = split('_', $j);
    @temp2 = split('\.', @geneReplicon[$n+1]);
    if(@temp[0] eq @temp2[0]) {
      do {
        $n++;
        say $out ("@geneNames[$n]\t@geneStart[$n]\t@geneEnd[$n]\t$geneLength\t0\t0\tNA");
        @temp = split('_', $j);
        @temp2 = split('\.', @geneReplicon[$n+1]);
      }until(@temp[0] ne @temp2[0])
    }
    close($in);
    close($out);
    system("sort -u -n -k2,2 $output > temp.txt");
    system("mv temp.txt $output");
    $content = '';
    open($in, '<', $output);
    while(<$in>) {
      $content = $content . $_;
    }
    close($in);
    open($out, '>', $output);
    say $out ("Gene_name\tStart_nt\tEnd_nt\tLength\tMotif_count\tNormalized_count\tExtent_methylation");
    print $out ($content);
    close($out);
  }
}
