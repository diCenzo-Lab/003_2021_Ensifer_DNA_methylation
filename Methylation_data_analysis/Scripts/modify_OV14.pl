#!usr/bin/perl
use 5.010;

while(<>) {
  if(/\#\#/) {
    print($_);
  }
  elsif(/m6A/) {
    chomp;
    if(/motif=/) {
      say($_);
    }
    else {
      @line = split("\t", $_);
      @line2 = split(';', @line[8]);
      @line2[0] =~ s/context=//;
      $motif = substr(@line2[0], 19, 5);
      if($motif eq 'GAGTC' || $motif eq 'GAATC' || $motif eq 'GATTC' || $motif eq 'GACTC') {
        say("@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t@line[7]\tcontext=@line2[0];@line2[1];motif=GANTC;@line2[2];@line2[3];id=GANTC;@line2[4];@line2[5];@line2[6];@line2[7]");
      }
      else {
        $motif = substr(@line2[0], 17, 8);
        if($motif eq 'AACCGATG' || $motif eq 'ATCCGATG' || $motif eq 'AGCCGATG' || $motif eq 'ACCCGATG' || $motif eq 'TACCGATG' || $motif eq 'TTCCGATG' || $motif eq 'TGCCGATG' || $motif eq 'TCCCGATG') {
          say("@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t@line[7]\tcontext=@line2[0];@line2[1];motif=WNCCGATG;@line2[2];@line2[3];id=WNCCGATG;@line2[4];@line2[5];@line2[6];@line2[7]");
        }
      }
    }
  }
}
