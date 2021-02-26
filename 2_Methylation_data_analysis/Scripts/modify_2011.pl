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
        $motif1 = substr(@line2[0], 17, 4);
        $motif2 = substr(@line2[0], 26, 3);
        if($motif1 eq 'CGCA' && $motif2 eq 'GTG') {
          say("@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t@line[7]\tcontext=@line2[0];@line2[1];motif=CGCANNNNNGTG;@line2[2];@line2[3];id=CGCANNNNNGTG/CACNNNNNTGCG;@line2[4];@line2[5];@line2[6];@line2[7]");
        }
        else {
          $motif1 = substr(@line2[0], 19, 3);
          $motif2 = substr(@line2[0], 27, 4);
          if($motif1 eq 'CAC' && $motif2 eq 'TGCG') {
            say("@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t@line[7]\tcontext=@line2[0];@line2[1];motif=CACNNNNNTGCG;@line2[2];@line2[3];id=CGCANNNNNGTG/CACNNNNNTGCG;@line2[4];@line2[5];@line2[6];@line2[7]");
          }
        }
      }
    }
  }
  elsif(/m4C/) {
    chomp;
    if(/motif=/) {
      say($_);
    }
    else {
      @line = split("\t", $_);
      @line2 = split(';', @line[8]);
      @line2[0] =~ s/context=//;
      $motif = substr(@line2[0], 16, 7);
      if($motif eq 'ACGCCTC' || $motif eq 'GCGCCTC' ) {
        say("@line[0]\t@line[1]\t@line[2]\t@line[3]\t@line[4]\t@line[5]\t@line[6]\t@line[7]\tcontext=@line2[0];@line2[1];motif=RCGCCTC;@line2[2];@line2[3];id=RCGCCTC;@line2[4];@line2[5];@line2[6];@line2[7]");
      }
    }
  }
}
