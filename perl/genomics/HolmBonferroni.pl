#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile.txt> <alpha>\n" unless @ARGV==2;
my ($infile,$alpha)=@ARGV;

my @A;
open(IN,$infile) || die "can't open $infile\n";
while(<IN>) {
  chomp;
  if(/\d/) { push @A,$_ }
}
close(IN);

@A=sort {$a<=>$b} @A;
my $n=@A;
for(my $i=0 ; $i<$n ; ++$i) {
  my $k=$i+1;
  my $a=$A[$i];
  if($a>$alpha/($n+1-$k)) {
    if($i==0) { print "no rejections possible\n"; exit }
    my $c=$A[$i-1];
    print "corrected alpha = $c\nreject $i tests\n";
    exit;
  }
}
print "reject all\n";

