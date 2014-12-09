#!/usr/bin/perl
#convert a GFF gene structure into a simple graph: 1 where there is an exon, 0 otherwise

use warnings;
use strict;

die "gff-to-graph.pl <gff>" unless @ARGV==1;

my $file = shift;
open(GFF, $file);
my @graph = ();
while(<GFF>)
{
  if (/^#/) {next;}
  my @lineSplit = split(/\t+/);
  my ($begin, $end) = ($lineSplit[3], $lineSplit[4]);
  for (my $i = $begin; $i <= $end; $i++)
  {
    $graph[$i] = 1;
  }
}

for (my $j = 0; $j < @graph; $j++)
{
  if (!(defined($graph[$j])))
  {
    $graph[$j] = 0;
  }
  print("$j $graph[$j]\n");
}
