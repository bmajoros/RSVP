#!/usr/bin/perl

use warnings;
use strict;

die "evaluate-exons.pl <prediction-gff> <annotation-gff>" unless @ARGV==2;

my ($pred, $ann) = @ARGV;
my ($numExonsTotal, $numExonsMatched) = (0,0);
open(ANN, "$ann") or die "unable to open annotation file";
my %exonHash = ();
while (<ANN>)
{
  if ($_ =~ /^#/) {next;}
  my @lineSplit = split(/\s+/, $_);
  my $exon = "$lineSplit[3] $lineSplit[4]";
  $exonHash{$exon} = 0;
}
close(ANN);
open(PRED, "$pred") or die "unable to open prediction file";
while (<PRED>)
{
  if ($_ =~ /^#/) {next;}
  my @lineSplit = split(/\s+/, $_);
  my $exon = "$lineSplit[3] $lineSplit[4]";

  my $annExonVal = $exonHash{$exon};
  if (defined($annExonVal))
  {
    $exonHash{$exon} = 1;
  }
}
close(PRED);
for my $annoExon (keys(%exonHash))
{
  $numExonsTotal++;
  if ($exonHash{$annoExon} == 1)
  {
    $numExonsMatched++;
  }
}
my $sn = $numExonsMatched/$numExonsTotal;
print("$sn\n");
#print("$numExonsMatched $numExonsTotal\n");
