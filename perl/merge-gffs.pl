#!/usr/bin/perl
#merge GFF predictions into one big file, and translate internal coords
#to real coords
#
#need output GFF directory and chunk directory with coords files

use warnings;
use strict;

die "merge-gffs.pl <gff-directory> <chunk-directory> [-comments]" unless @ARGV == 2 || @ARGV == 3;

my ($gffDir, $chunkDir) = ($ARGV[0], $ARGV[1]);
my $comments;
if (exists($ARGV[2]) && $ARGV[2] eq "-comments")
{
  $comments = 1;
}
else
{
  $comments = 0;
}
my $gffPrint = `ls $gffDir/*.gff`;
my @gffPrintSplit = split(/\s+/, $gffPrint);

my @gffList = ();
for my $gffName (@gffPrintSplit)
{
  if (!($gffName =~ /forward/ || $gffName =~ /reverse/))
  {
    push(@gffList, $gffName);
  }
}

for my $gff (@gffList)
{
  my $chunkNum;
  if ($gff =~ /(\d+).gff/)
  {
    $chunkNum = $1;
  }
  else
  {
    next;
  }
  my $coordsFile = "$chunkDir/$chunkNum.coords.gff";
  open(GFF, $gff) or die "Can't find gff prediction for chunk $chunkNum";
  open(COORD, $coordsFile) or die "Can't find coords.gff file for chunk $chunkNum: file is $coordsFile";

  my $coordLine = <COORD>;
  chomp($coordLine);
  my @coordLineSplit = split(/\s+/, $coordLine);
  my $chromosome = $coordLineSplit[0];
  my $chunkStart = $coordLineSplit[3];
  my $coordTranscriptString = $coordLineSplit[8];
  my $gene;
  if ($coordTranscriptString =~ /transcript_id=(\w+).\d+/)
  {
    $gene = $1;
  }
  elsif ($coordTranscriptString =~ /gene_id=(\w+)/)
  {
    $gene = $1;
  }
  else
  {
    die "could not find gene from coords.gff file for chunk $chunkNum";
  }

  while(<GFF>)
  {
    my $currentLine = $_;
    if ($currentLine =~ /^#/)
    {
      if ($comments == 1)
      {
        print($currentLine);
      }
      next;
    }
    my @currentLineSplit = split(/\s+/, $currentLine);
    my $feature;
    if ($currentLineSplit[2] =~ /exon/)
    {
      $feature = "exon";
    }
    elsif ($currentLineSplit[2] =~ /intron/)
    {
      $feature = "intron";
    }
    else
    {
      die "found feature that is not exon or intron: $feature";
    }
    my ($featureStart, $featureEnd, $gffTranscriptString) = ($currentLineSplit[3], $currentLineSplit[4], $currentLineSplit[8]);
    my $gffSupportString = ($currentLineSplit[9]);
    my $newStart = $featureStart + $chunkStart - 1; #because both coordinates are 1-based
    my $newEnd = $featureEnd + $chunkStart - 1;
    my ($transcriptNum, $parseNum);
    if ($gffTranscriptString =~ /transcript_id=(\d+)/)
    {
      $transcriptNum = $1;
    }
    else
    {
      die "transcript_id error";
    }
    if ($gffTranscriptString =~ /parse=(\d+)/)
    {
      $parseNum = $1;
    }
    my $newTranscriptLabel;
    if (defined($parseNum))
    {
      $newTranscriptLabel = "transcript_id \"$gene.$transcriptNum.$parseNum\"; gene_id \"$gene\";";
    }
    else
    {
      $newTranscriptLabel = "transcript_id \"$gene.$transcriptNum\"; gene_id \"$gene\";";
    }
    my $newSupportLabel;
    if (defined($gffSupportString))
    {
      my ($numSupport, $numTotalSupport);
      if ($gffSupportString =~ /support=([\d-]+);fullSupport=([\d-]+)/)
      {
        ($numSupport, $numTotalSupport) = ($1, $2);
      }
      else
      {
        die "support string error $gene";
      }
      $newSupportLabel = " support \"$numSupport\"; fullSupport \"$numTotalSupport\";";
    }
    else
    {
      $newSupportLabel = "";
    }
    my $newLine = "$chromosome\t$currentLineSplit[1]\t$feature\t$newStart\t$newEnd\t$currentLineSplit[5]\t$currentLineSplit[6]\t$currentLineSplit[7]\t$newTranscriptLabel$newSupportLabel\n";
    print($newLine);
  }

  close(GFF);
  close(COORD);
}
