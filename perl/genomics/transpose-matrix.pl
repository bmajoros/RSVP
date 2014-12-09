#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile>\n" unless @ARGV==1;
my ($infile)=@ARGV;

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
my $m=@{$M[0]};
for(my $i=0 ; $i<$m ; ++$i) {
  for(my $j=0 ; $j<$n ; ++$j) {
    my $x=$M[$j]->[$i];
    print "$x";
    if($j<$n-1) {print "\t"}
  }
  print "\n";
}
