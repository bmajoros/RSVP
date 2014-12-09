#!/usr/bin/perl
#Make a pdiff file specified by Song's description

use warnings;
use strict;

die "make-pdiff.pl <annotation-gff> <prediction-gff>" unless @ARGV==2;

my ($ann, $pred) = @ARGV;

my ($annStart, $annEnd) = (9**9**9, -9**9**9);
my ($predStart, $predEnd) = (9**9**9, -9**9**9);

my @annStructure = ();
my @predStructure = ();

#first pass: get start and end to determine offset

open(ANN1, $ann);

while(<ANN1>)
{
  if (/^#/) { next; }
  my @lineSplit = split(/\s+/);
  my ($fStart, $fEnd) = ($lineSplit[3], $lineSplit[4]);
  if ($fStart < $annStart) { $annStart = $fStart; }
  if ($fEnd > $annEnd) { $annEnd = $fEnd; }
}

close(ANN1);

open(PRED1, $pred);

while(<PRED1>)
{
  if (/^#/) { next; }
  my @lineSplit = split(/\s+/);
  my ($fStart, $fEnd) = ($lineSplit[3], $lineSplit[4]);
  if ($fStart < $predStart) { $predStart = $fStart; }
  if ($fEnd > $predEnd) { $predEnd = $fEnd; }
}

close(PRED1);

my $offset = ($predStart < $annStart ? $predStart : $annStart);

open(ANN, $ann);

while(<ANN>)
{
  if (/^#/) { next; }
  my @lineSplit = split(/\s+/);
  my ($fStart, $fEnd) = ($lineSplit[3], $lineSplit[4]);
  for (my $i = $fStart; $i <= $fEnd; $i++)
  {
    $annStructure[$i - $offset] = 1;
  }
}

close(ANN);

for (my $j = 0; $j < @annStructure; $j++)
{
  if (!defined($annStructure[$j]))
  {
    $annStructure[$j] = 0;
  }
}

open(PRED, $pred);

while(<PRED>)
{
  if (/^#/) { next; }
  my @lineSplit = split(/\s+/);
  my ($fStart, $fEnd) = ($lineSplit[3], $lineSplit[4]);
  for (my $k = $fStart; $k <= $fEnd; $k++)
  {
    $predStructure[$k - $offset] = 1;
  }
}

close(PRED);

for (my $l = 0; $l < @predStructure; $l++)
{
  if (!defined($predStructure[$l]))
  {
    $predStructure[$l] = 0;
  }
}
my $predLength = @predStructure;
my $annLength = @annStructure;
my $maxLength = ($predLength > $annLength ? $predLength : $annLength);

my ($predPhase, $annPhase) = (0,0);

for (my $m = 0; $m < $maxLength; $m++)
{
  my $pdiffState;
  my $predRawState = $predStructure[$m];
  my $annRawState = $annStructure[$m];
  my ($predState, $annState);
  my $coord = $m + $offset;
  if ($coord < $predStart || $coord > $predEnd) { $predState = "UTR"; }
  elsif ($predRawState == 0) { $predState = "intron"; }
  elsif ($predRawState == 1) { $predState = "exon"; }
  if ($coord < $annStart || $coord > $annEnd) { $annState = "UTR"; }
  elsif ($annRawState == 0) { $annState = "intron"; }
  elsif ($annRawState == 1) { $annState = "exon"; }

  if ($annState eq "exon" && $predState eq "exon")
  {
    if ($annPhase == $predPhase)
    {
      $pdiffState = "CCI";
    }
    else
    {
      $pdiffState = "CCO";
    }
  }
  if ($annState eq "intron" && $predState eq "exon")
  {
    $pdiffState = "IC";
  }
  if ($annState eq "exon" && $predState eq "intron")
  {
    $pdiffState = "CI";
  }
  if ($annState eq "intron" && $predState eq "intron")
  {
    $pdiffState = "I";
  }
  if ($annState eq "exon" && $predState eq "UTR")
  {
    $pdiffState = "CU";
  }
  if ($annState eq "UTR" && $predState eq "exon")
  {
    $pdiffState = "UC";
  }
  if ($annState eq "intron" && $predState eq "UTR")
  {
    $pdiffState = "IU";
  }
  if ($annState eq "UTR" && $predState eq "intron")
  {
    $pdiffState = "UI";
  }
  if ($annState eq "UTR" && $predState eq "UTR")
  {
    $pdiffState = "U";
  }

  if ($annState eq "exon")
  {
    if ($annPhase == 2)
    {
      $annPhase = 0;
    }
    else
    {
      $annPhase++;
    }
  }
  if ($predState eq "exon")
  {
    if ($predPhase == 2)
    {
      $predPhase = 0;
    }
    else
    {
      $predPhase++;
    }
  }
  
  print("$coord $pdiffState\n");
}
