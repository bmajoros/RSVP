#!/usr/bin/perl
#only works on single-parse, single-gene outputs right now

use warnings;
use strict;
use Switch;
use Translation;

die "translate-protein.pl <*.gff> <*.fasta>" unless @ARGV==2;

my ($gff, $fasta) = @ARGV;

open(GFF, $gff);
my @gff = <GFF>;
close(GFF);

open(FASTA, $fasta);
my @fasta = <FASTA>;
close(FASTA);

my $sequence;
my @structure;
my $strand;

#make list where at each index, there is a 1 if that coordinate has an exon base
#this code was lifted from compare-predictions-asta-new.pl
for (my $i = 0; $i < @gff; $i++)
{
  my $currentFeature = $gff[$i];
  if ($currentFeature =~ /^#/) {next};
  my @featureSplit = split(/\s+/,$currentFeature);
  my ($begin, $end) = ($featureSplit[3], $featureSplit[4]);
  my $currentStrand = $featureSplit[6];
  if (!defined($strand)) { $strand = $currentStrand; }
  if (!($strand eq $currentStrand)) { die "error: all features are not on the same strand"; }
  for (my $j = $begin; $j <= $end; $j++)
  {
    $structure[$j] = 1;
  }
}
for (my $k = 0; $k < @structure; $k++) #set all of the other elements to 0
{
  if (!defined($structure[$k]))
  {
    $structure[$k] = 0;
  }
}

for (my $f = 1; $f < @fasta; $f++)
{
  my $currentLine = $fasta[$f];
  chomp($currentLine);
  $sequence = "${sequence}${currentLine}";
}

my @fwdProtein = ();
my $fwdPhase = 0;
my @fwdBase = ();
my @fwdCoord = ();

for (my $i = 1; $i < @structure; $i++)
{
  if ($structure[$i] == 0) {next;}
  my $curBase = substr($sequence, $i-1, 1); #structure is 1-based but fasta string is 0-based
  $fwdBase[$fwdPhase] = $curBase;
  $fwdCoord[$fwdPhase] = $i;
  if ($fwdPhase == 2)
  {
    #translate 3 bases to amino acid here
    my $baseString = "$fwdBase[0]$fwdBase[1]$fwdBase[2]";
    my $bsRef = \$baseString;
    my $aminoAcid = Translation::translate($bsRef);
    my $aaString = "$aminoAcid\t$fwdCoord[0]\t$fwdCoord[1]\t$fwdCoord[2]";
    push(@fwdProtein, $aaString);
    $fwdPhase = 0;
  }
  else
  {
    $fwdPhase++;
  }
}

my @revProtein = ();
my $revPhase = 0;
my @revCompBase = ();
my @revCoord = ();

for (my $j = @structure - 1; $j > 0; $j--)
{
  if ($structure[$j] == 0) {next;}
  my $curBase = substr($sequence, $j-1, 1);
  my $compBase = comp($curBase);
  $revCompBase[$revPhase] = $compBase;
  $revCoord[$revPhase] = $j;
  if ($revPhase == 2)
  {
    #translate 3 bases to amino acid here
    my $baseString = "$revCompBase[0]$revCompBase[1]$revCompBase[2]";
    my $bsRef = \$baseString;
    my $aminoAcid = Translation::translate($bsRef);
    my $aaString = "$aminoAcid\t$revCoord[0]\t$revCoord[1]\t$revCoord[2]";
    push(@revProtein, $aaString);
    $revPhase = 0;
  }
  else
  {
    $revPhase++;
  }
}

my @curProtein;
if ($strand eq "+")
{
  @curProtein = @fwdProtein;
}
elsif ($strand eq "-")
{
  @curProtein = @revProtein;
}

print("#Protein strand: $strand\n");

for (my $x = 0; $x < @curProtein; $x++)
{
  print("$curProtein[$x]\n");
}

sub comp
{
  my $oldBase = shift;
  my $compBase;
  switch($oldBase)
  {
    case "A" { $compBase = "T" }
    case "T" { $compBase = "A" }
    case "C" { $compBase = "G" }
    case "G" { $compBase = "C" }
  }
  return $compBase;
}
