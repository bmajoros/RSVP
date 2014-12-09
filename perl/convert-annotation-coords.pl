#!/usr/bin/perl
#translate internal coords of annotations to real coords
#based on convert-gff-coords.pl
#need chunk directory with coords files

use warnings;
use strict;

die "convert-annotation-coords.pl <chunk-directory> <out-directory>" unless @ARGV == 2;

my ($chunkDir, $outDir) = ($ARGV[0], $ARGV[1]);
my $gffPrint = `ls $chunkDir/*.gff`;
my @gffPrintSplit = split(/\s+/, $gffPrint);

my @gffList = ();
for my $gffName (@gffPrintSplit)
{
  if (!($gffName =~ /coords/))
  {
    push(@gffList, $gffName);
  }
}

for my $gff (@gffList)
{
  my $gffFile;
  if ($gff =~ /\/(\w+\.gff)/) {$gffFile = $1;}
  open(OUT, "> $outDir/$gffFile") or die "error--could not create output file";
  my ($chunkNum) = ($gff =~ /(\d+).gff/);
  my $coordsFile = "$chunkDir/$chunkNum.coords.gff";
  open(GFF, $gff) or die "Can't find gff annotation for chunk $chunkNum";
  open(COORD, $coordsFile) or die "Can't find coords.gff file for chunk $chunkNum: file is $coordsFile";

  my $coordLine = <COORD>;
  chomp($coordLine);
  my @coordLineSplit = split(/\s+/, $coordLine);
  my $chromosome = $coordLineSplit[0];
  my $chunkStart = $coordLineSplit[3];

  while(<GFF>)
  {
    my $currentLine = $_;
    chomp($currentLine);
    if ($currentLine =~ /^#/)
    {
      next;
    }
    my @currentLineSplit = split(/\t+/, $currentLine);
    if (!($currentLineSplit[2] =~ /exon/))
    {
      die("non-exon feature found");
    }
    my ($featureStart, $featureEnd) = ($currentLineSplit[3], $currentLineSplit[4]);
    my $newStart = $featureStart + $chunkStart - 1; #because both coordinates are 1-based
    my $newEnd = $featureEnd + $chunkStart - 1;
    my $newLine = "$chromosome\t$currentLineSplit[1]\texon\t$newStart\t$newEnd\t$currentLineSplit[5]\t$currentLineSplit[6]\t$currentLineSplit[7]\t$currentLineSplit[8]\n";
    print(OUT $newLine);
  }
  close(OUT);
  close(GFF);
  close(COORD);
}
