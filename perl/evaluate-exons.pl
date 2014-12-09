#!/usr/bin/perl

use warnings;
use strict;

die "evaluate-exons.pl <out-dir> <data-dir> <extension>\n" unless @ARGV==3;
my ($out, $data, $ext) = @ARGV;

my $gffs = `ls $out/*.$ext`;
my @gffList = split(/\s+/,$gffs);
my ($numExonsTotal, $numExonsMatched) = (0,0);
my $numPreds=0;
for my $gff (@gffList)
{
  my $gffNum;
  if ($gff =~ /(\d+).$ext$/)
  {
    $gffNum = $1;
  }
  else
  {
    next;
  }
  open(ANN, "$data/$gffNum.gff") or die "unable to open annotation $gffNum";
  my %exonHash = ();
  while (<ANN>)
  {
    if ($_ =~ /^#/) {next;}
    my @lineSplit = split(/\s+/, $_);
    next unless @lineSplit>8;

    next unless $lineSplit[2] eq "internal-exon"; ###

    my $exon = "$lineSplit[3] $lineSplit[4]";
    $exonHash{$exon} = 0;
  }
  close(ANN);
  open(PRED, "$out/$gffNum.$ext") or die "unable to open prediction $gffNum";
  while (<PRED>)
  {
    if ($_ =~ /^#/) {next;}
    my @lineSplit = split(/\s+/, $_);
    next unless @lineSplit>8;

    my $exon = "$lineSplit[3] $lineSplit[4]";
    my $annExonVal = $exonHash{$exon};
    if (defined($annExonVal))
    {
      $exonHash{$exon} = 1;
    }
    ++$numPreds;
  }
  for my $annoExon (keys(%exonHash))
  {
    $numExonsTotal++;
    if ($exonHash{$annoExon} == 1)
    {
      $numExonsMatched++;
    }
  }
  close(PRED);
}
my $sn = $numExonsMatched/$numExonsTotal;
print("$numExonsMatched / $numExonsTotal = $sn ($numPreds predictions)\n");
#print("$numExonsMatched $numExonsTotal\n");
