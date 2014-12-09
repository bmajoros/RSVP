#!/usr/bin/perl
#Delete predictions with a parse score of -Infinity from the output GFF.

use warnings;
use strict;

die "remove-inf.pl <gff>" unless @ARGV==1;

my $gff = shift;
my @infParses = ();
open(FIRST_TIME, $gff);
while(<FIRST_TIME>)
{
  if ($_ =~ /# parse (\d+) score: -Infinity/)
  {
    push(@infParses, $1);
  }
}
close(FIRST_TIME);
open(SECOND_TIME, $gff);
while(<SECOND_TIME>)
{
  my $toPrint = 1;
  for (my $i = 0; $i < @infParses; $i++)
  {
    my $parseNum = $infParses[$i];
    if ($_ =~ /parse=$parseNum;/ || $_ =~ /parse $parseNum score:/)
    {
      $toPrint = 0;
    }
  }
  if ($toPrint == 1)
  {
    print($_);
  }
}
