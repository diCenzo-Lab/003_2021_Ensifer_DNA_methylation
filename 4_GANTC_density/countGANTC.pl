#!usr/bin/perl
use 5.010;

$count = 0;
$GANTC = 0;
$GC = 0;
$nt = 0;
$temp_file = 'temporary.txt';
open($out, '>', $temp_file);
while(<>) {
  chomp;
  if(/>/) {
    $count++;
    print $out ("\n") unless($count == 1);
  }
  else {
    print $out ("$_");
  }
}
close($out);

open($in, '<', $temp_file);
while(<$in>) {
  @line = split('', $_);
  $length = scalar(@line);
  for($n=0; $n < $length - 4; $n++) {
    $test = @line[$n] . @line[$n+1] . @line[$n+2] . @line[$n+3] . @line[$n+4];
    if($test eq 'GAATC' || $test eq 'GAGTC' || $test eq 'GATTC' || $test eq 'GACTC') {
      $GANTC++;
    }
  }
  for($n=0; $n < $length; $n++) {
    $nt++;
    if(@line[$n] eq 'G' || @line[$n] eq 'C') {
      $GC++;
    }
  }
}
close($in);
$density = 1000 * $GANTC / $nt;
$GC_content = 100 * $GC / $nt;
say("\t$GANTC\t$nt\t$density\t$GC_content");
unlink($temp_file);
