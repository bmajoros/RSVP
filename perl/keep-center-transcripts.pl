#!/usr/bin/perl
#Reads in a GFF prediction file from Decoder and deletes all but the center
#transcript (or the closest transcript left of center) in each parse
use strict;

die "keep-center-transcripts.pl <decoder-gff>" unless @ARGV==1;

my $decoderFile = shift(@ARGV);
open(DECODER_FILE, $decoderFile);
my @decoderGff = <DECODER_FILE>;
my $outputGff = keepCenterTranscripts(\@decoderGff);
print(@$outputGff);

sub keepCenterTranscripts
{
  my $gffref = shift;
  my @gff = @$gffref;
  my @outputgff = ();
  #first pass: count number of parses, get min begin and max end
  my ($numParses, $currentParse, $minBegin, $maxEnd) = (0,0,9**9**9,-9**9**9);
  for (my $i = 0; $i < @gff; $i++)
  {
    if ($gff[$i] =~ /^#/)
    {
      next;
    }
    my @splitLine = split(/\s+/, $gff[$i]);
    my ($begin, $end) = ($splitLine[3], $splitLine[4]);
    my $parseNum = 0;
    if ($gff[$i] =~ /transcript_id \"*\w+\.\d+\.(\d+)\"*;/)
    {
      $parseNum = $1;
    }
    else
    {
      die("parseNum error in prediction file");
    }
    if ($parseNum != $currentParse)
    {
      $numParses++;
      $currentParse = $parseNum;
    }
    if ($begin < $minBegin)
    {
      $minBegin = $begin;
    }
    if ($end > $maxEnd)
    {
      $maxEnd = $end;
    }
  }
  my $center = ($maxEnd-$minBegin)/2 + $minBegin;
  #second pass: find the center (relative to minBegin and maxEnd) transcript in each parse
  #find the beginning and end of each transcript in the current parse, and pick the right one
  my @centerTranscripts;
  my $currentParseB = 0;
  my @tBegin = (); #keeps track of the beginning of each transcript in the current parse
  my @tEnd = (); #keeps track of the end of each transcript in the current parse
  for (my $j = 0; $j < @gff; $j++) #@gff; $j++)
  {
    next if $gff[$j] =~ /^#/;
    my @splitLine = split(/\s+/, $gff[$j]);
    my ($begin, $end, $transcriptId) = ($splitLine[3], $splitLine[4], $splitLine[8]);
    my ($transcriptNum, $parseNum);
    if ($gff[$j] =~ /transcript_id \"*\w+\.(\d+)\.(\d+)\"*;/)
    {
      ($transcriptNum, $parseNum) = ($1, $2);
    }
    else
    {
      print("transcript and parseNum error in prediction file\n");
    }
    if ($parseNum !=  $currentParseB)
    {
      $currentParseB = $parseNum;
      @tBegin = ();
      @tEnd = ();
    }
    if (!defined($tBegin[$transcriptNum]) || $tBegin[$transcriptNum] > $begin)
    {
      $tBegin[$transcriptNum] = $begin;
    }
    if (!defined($tEnd[$transcriptNum]) || $tEnd[$transcriptNum] < $end)
    {
      $tEnd[$transcriptNum] = $end;
    }
    #pick the center parse or the closest parse left of center
    for (my $a = 1; $a < @tBegin; $a++)
    {
      if ($tBegin[$a] <= $center && $tEnd[$a] >= $center)
      {
        $centerTranscripts[$currentParseB] = $a;
      }
      elsif ($tEnd[$a] <= $center && ($a + 1 >= @tBegin || $tBegin[$a + 1] > $center))
      {
        $centerTranscripts[$currentParseB] = $a;
      }
    }
  }
  #third pass: delete everything but the center transcript in each parse
  for (my $k = 0; $k < @gff; $k++)
  {
    if ($gff[$k] =~ /^#/)
    {
      next;
    }
    my @splitLine = split(/\s+/, $gff[$k]);
    my ($begin, $end, $transcriptId) = ($splitLine[3], $splitLine[4], $splitLine[8]);
    my ($transcriptNum, $parseNum);
    if ($gff[$k] =~ /transcript_id \"*\w+\.(\d+)\.(\d+)\"*;/)
    {
      ($transcriptNum, $parseNum) = ($1, $2);
    }
    else
    {
      print("transcript and parseNum error in prediction file\n");
    }
    my $centerNum = $centerTranscripts[$parseNum];
    if ($centerNum == $transcriptNum)
    {
      push(@outputgff, $gff[$k]);
    }
  }
  return \@outputgff;
}

