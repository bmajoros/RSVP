#!/usr/bin/perl
#Print out the top parse scores in a tab delimited format:
#chunk1.parse1	chunk1.parse2	...
#chunk2.parse1	chunk2.parse2	...

use warnings;
use strict;

die "get-score-dist.pl <out-dir>" unless @ARGV==1;

my ($outDir) = @ARGV;
my @files = `ls $outDir/*.gff`;
for my $file (@files)
{
  open(GFF, "$file");
  my @gff = <GFF>;
  close(GFF);
  my @scores = ();
  for my $line (@gff)
  {
    next unless $line =~ /# parse \d+ score: (-\d+\.\d+)/;
    push(@scores, "$1");
  }
  my $firstScore = $scores[0];
  for (my $i=1; $i<@scores; $i++)
  {
    my $score = $scores[$i];
    my $diffScore = $score-$scores[$i-1];
    my $normScore = $diffScore/$firstScore;
    print "$normScore\t";
  }
  print "\n";
}
