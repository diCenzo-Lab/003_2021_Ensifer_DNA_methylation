#!usr/bin/perl
use 5.010;

$n = 0;
while(<>) {
  $n++;
  chomp;
  @line = split("\t", $_);
  foreach $i (@line) {
    if($i ne '') {
      print("$i\t$n\n");
    }
  }
}
