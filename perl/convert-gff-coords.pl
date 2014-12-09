#!/usr/bin/perl
#Translates internal coords of GFF predictions to real coords and
#formats attribute fields according to the GTF standard.

#need prediction GFF directory and chunk directory with coords files
#converts in place

use warnings;
use strict;

die "convert-gff-coords.pl <gff-directory> <ext> <chunk-directory> <out-directory> [-comments]" unless @ARGV == 4 || @ARGV == 5;

my ($gffDir, $ext, $chunkDir, $outDir) = @ARGV;
my $comments;
if (exists($ARGV[3]) && $ARGV[3] eq "-comments")
{
  $comments = 1;
}
else
{
  $comments = 0;
}
my $gffPrint = `ls $gffDir/*.$ext`;
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
  my $gffFile;
  if ($gff =~ /\/(\w+\.$ext)/) {$gffFile = $1;}
  open(OUT, "> $outDir/$gffFile") or die "error--could not create output file";
  my ($chunkNum) = ($gff =~ /(\d+).$ext/);
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
    chomp($currentLine);
    if ($currentLine =~ /^#/)
    {
      if ($comments == 1)
      {
        print($currentLine);
      }
      next;
    }
    my @currentLineSplit = split(/\s+/, $currentLine);
    next unless @currentLineSplit>8;
    my ($featureStart, $featureEnd) = ($currentLineSplit[3], $currentLineSplit[4]);
    my $newStart = $featureStart + $chunkStart - 1; #because both coordinates are 1-based
    my $newEnd = $featureEnd + $chunkStart - 1;
    my ($transcriptNum, $parseNum);
    my ($gffTranscriptString, $gffSupportString) = ($currentLineSplit[8], $currentLineSplit[9]);
    if ($gffTranscriptString =~ /transcript_id="([^";]+)/)
    {
      $transcriptNum = $1;
    }
    elsif ($gffTranscriptString =~ /transcript_id=([^";]+)/)
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
      die "unknown feature type: $currentLineSplit[2]";
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
    print(OUT $newLine);
  }
  close(OUT);
  close(GFF);
  close(COORD);
}
