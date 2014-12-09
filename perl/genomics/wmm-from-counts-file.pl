#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile> <pseudocount>\n" unless @ARGV==2;
my ($infile,$pseudo)=@ARGV;

my @M;
open(IN,$infile) || die "can't open file: $infile\n";
while(<IN>) {
  chomp;
  if(/(\S.*\S)/) {
    my @fields=split/\s+/,$1;
    my $f=[];
    @$f=@fields;
    push @M,$f;
  }
}
close(IN);
my $n=@M;
print "WMM\nPROMOTER\n";
print "0 $n 4\n$n 0 0 +\n";
for(my $i=0 ; $i<$n ; ++$i) {
  my $row=$M[$i];
  my $sum=0;
  for(my $j=0 ; $j<4 ; ++$j) {
    $row->[$j]+=$pseudo;
    $sum+=$row->[$j];
  }
  for(my $j=0 ; $j<4 ; ++$j) {
    my $x=log($row->[$j]/$sum);
    print "$x\t";
  }
  print "\n";
}

