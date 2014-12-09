#!/usr/bin/perl
#Compare predictions to GFF annotations using AStalavista

use warnings;
use strict;

sub format;
sub compare;
sub splitIntoParses;
sub getNucleotides;

die "compare-predictions-asta.pl <annotation-gff> <prediction-gff> <astalavista-path>\n" unless @ARGV==3;

#print("#Note: prediction and annotation numbers are 1-based and correspond to their positions in
##their respective GFF files. For instance, the top prediction is prediction 1, the second
##is prediction 2, etc.\n");

my ($aFile, $pFile, $astapath) = @ARGV;

open(A_FILE, $aFile);
open(P_FILE, $pFile);

#Arrays of strings (GFF lines) for prediction and annotation
my @annotation = (<A_FILE>);
my @prediction = (<P_FILE>);

my $aSplitRef = splitIntoParses(\@annotation);
my $pSplitRef = splitIntoParses(\@prediction);

for (my $i = 0; $i < @$pSplitRef; $i++)
{
  my $currentPred = $pSplitRef->[$i];
  #figure out how far to compare
  my $pFirstLine = $currentPred->[0]; #either the first or last line will have the ending feature
  if (!(defined($pFirstLine))) {die "error: prediction has no lines. Corresponding annotation file: $aFile"}
  my @pFirstLineSplit = split(/\s+/,$pFirstLine);
  my ($pFirstLineBegin, $pFirstLineEnd) = ($pFirstLineSplit[3], $pFirstLineSplit[4]);
  my $pLastLine = $currentPred->[@$currentPred-1];
  my @pLastLineSplit = split(/\s+/,$pLastLine);
  my ($pLastLineBegin, $pLastLineEnd) = ($pLastLineSplit[3], $pLastLineSplit[4]);
  my $endBound = ($pFirstLineEnd > $pLastLineEnd ? $pFirstLineEnd : $pLastLineEnd); #only compare to end of prediction
  #compare pairwise to all the annotations
  my $minDiff = 9**9**9; #infinity
  my $minDiffAnn;
  my $minDiffAnnNum;
  for (my $j = 0; $j < @$aSplitRef; $j++)
  {
    my $currentAnn = $aSplitRef->[$j];
    #find starting coordinate for getNucleotides() arrays
    my $aFirstLine = $currentAnn->[0]; #same as above
    my @aFirstLineSplit = split(/\s+/,$aFirstLine);
    my ($aFirstLineBegin, $aFirstLineEnd) = ($aFirstLineSplit[3], $aFirstLineSplit[4]);
    my $aLastLine = $currentAnn->[@$currentAnn-1];
    my @aLastLineSplit = split(/\s+/,$aLastLine);
    my ($aLastLineBegin, $aLastLineEnd) = ($aLastLineSplit[3], $aLastLineSplit[4]);
    my $pBegin = ($pFirstLineBegin < $pLastLineBegin ? $pFirstLineBegin : $pLastLineBegin);
    my $aBegin = ($aFirstLineBegin < $aLastLineBegin ? $aFirstLineBegin : $aLastLineBegin);
    my $offset = ($pBegin < $aBegin ? $pBegin : $aBegin);

    my $predNucleotides = getNucleotides($currentPred, $offset);
    my $annNucleotides = getNucleotides($currentAnn, $offset);
    my $differences = 0;
    for (my $k = 0; $k <= $endBound - $offset; $k++)
    {
      if ($k >= @$predNucleotides || $k >= @$annNucleotides)
      {
        $differences += 1;
      }
      elsif ($predNucleotides->[$k] != $annNucleotides->[$k])
      {
        $differences += 1;
      }
    }
    if ($differences < $minDiff)
    {
      $minDiff = $differences;
      $minDiffAnn = $currentAnn;
      $minDiffAnnNum = $j;
    }
  }
  #get the transcript IDs of the current annotation and prediction
  my ($currentPredID, $currentAnnID);
  if ($pFirstLine =~ /transcript_id[\s=]\"*(\w*\.*\d+\.*\d*)\"*;/)
  {
    $currentPredID = $1;
  }
  else
  {
    die("Prediction transcript ID formatted incorrectly\n");
  }
  if ($minDiffAnn->[0] =~ /transcript_id[\s=]\"*(\w*\.*\d+)\"*;/)
  {
    $currentAnnID = $1;
  }
  else
  {
    die("Annotation transcript ID formatted incorrectly\n");
  }
  #compare the prediction with the closest annotation
  my $afRef = &format($minDiffAnn, "Temp", 1);
  my $pfRef = &format($currentPred, "Temp", 2);
  my ($predNum, $annNum) = ($i + 1, $minDiffAnnNum + 1);
  #print("Comparing prediction $predNum to annotation $annNum: ");
  print("$currentPredID\t$currentAnnID\t");
  compareAsta($afRef, $pfRef);
  print("\t");
  compareEnds($afRef, $pfRef);
  print("\n");
}

#use AStalavista to find the alternative splice events in two GFF files
sub compareAsta
{
  my ($aRef, $pRef) = @_;
  open(TEMPGFF, "> compare-predictions-temp-1.gff");
  print(TEMPGFF @$aRef);
  print(TEMPGFF @$pRef);
  close(TEMPGFF);
  system("$astapath/bin/astalavista.sh -c asta -i compare-predictions-temp-1.gff -o compare-predictions-temp-2.gff +ext 2>/dev/null");
  system("rm compare-predictions-temp-1.gff");
  system("gunzip compare-predictions-temp-2.gff.gz");
  open(ASTAFILE, "compare-predictions-temp-2.gff");
  my @astafile = <ASTAFILE>;
  close(ASTAFILE);
  system("rm compare-predictions-temp-2.gff");
  my $result = parseASTA(\@astafile);
  my @structureList = @$result;
  my @asEventList = ();
  if (@structureList == 0)
  {
    print("no difference ");
  }
  else
  {
    print("@structureList ");
  }
}

#Compare the first and last features in the two files to see if there are
#alternative start/stop codons
sub compareEnds
{
  my ($aRef, $pRef) = @_;
  my $afl = $aRef->[0];
  my $pfl = $pRef->[0];
  my $all = $aRef->[@$aRef-1];
  my $pll = $pRef->[@$pRef-1];
  my @aflSplit = split(/\s+/, $afl);
  my ($afb, $afe) = ($aflSplit[3], $aflSplit[4]);
  my @pflSplit = split(/\s+/, $pfl);
  my ($pfb, $pfe) = ($pflSplit[3], $pflSplit[4]);
  my @allSplit = split(/\s+/, $all);
  my ($alb, $ale) = ($allSplit[3], $allSplit[4]);
  my @pllSplit = split(/\s+/, $pll);
  my ($plb, $ple) = ($pllSplit[3], $pllSplit[4]);
  my $strand = $aflSplit[6];
  my ($flag1, $flag2) = (0,0);
  if (!($strand =~ /[\+\-]/))
  {
    die("error: improper strand info\n");
  }
  if ($afb != $pfb)
  {
    $flag1 = 1;
    if ($strand eq "+")
    {
      print("alt-start ");
    }
    elsif ($strand eq "-")
    {
      print("alt-stop ");
    }
  }
  if ($ale != $ple)
  {
    $flag2 = 1;
    if ($strand eq "+")
    {
      print("alt-stop ");
    }
    elsif ($strand eq "-")
    {
      print("alt-start ");
    }
  }
  if ($flag1 == 0 && $flag2 == 0)
  {
    print("no difference");
  }
}

#Parse the AStalavista output GFF file and get the alternative splice
#events in their encoding
sub parseASTA
{
  my $astaref = shift;
  my @astafile = @$astaref;
  my @structureList;
  for (my $i = 0; $i < @astafile; $i++)
  {
    my $structure;
    if ($astafile[$i] =~ /structure "(\S+)"/)
    {
      $structure = $1;
    }
    else
    {
      print("Error in AStalavista output pattern match\n");
    }
    push(@structureList, $structure);
  }
  return \@structureList;
}

##make an array for each parse/annotation where a[i] is 1 if i is in
#an exon, 0 otherwise
sub getNucleotides
{
  my ($tRef, $offset) = @_;
  my ($currentBegin, $currentEnd);
  my @nucleotideList = ();

  for (my $i = 0; $i < @$tRef; $i++) #set all of the elements in an exon to 1
  {
    my $currentFeature = $tRef->[$i];
    my @featureSplit = split(/\s+/,$currentFeature);
    my ($begin, $end) = ($featureSplit[3], $featureSplit[4]);
    if (!(defined($offset))) {die "error: offset not defined (feature $currentFeature)"}
    if ($begin < $offset) {die "error in getNucleotides: offset too high"}
    for (my $j = $begin; $j <= $end; $j++)
    {
      $nucleotideList[$j - $offset] = 1;
    }
  }

  for (my $k = 0; $k < @nucleotideList; $k++) #set all of the other elements to 0
  {
    if (!defined($nucleotideList[$k]))
    {
      $nucleotideList[$k] = 0;
    }
  }
  return \@nucleotideList;
}

#return an array of parses, each itself an array of lines
#need to change this if we're considering parses with multiple transcripts
sub splitIntoParses
{
  my $gffRef = shift;
  my @gff = @$gffRef;
  my @parses = ();
  my @currentParse = ();
  my $currentTranscript = "-1";
  for (my $i = 0; $i < @gff; $i++)
  {
    my $currentLine = $gff[$i];
    my $transcriptId;
    if ($currentLine =~ /^#/)
    {
      next;
    }
    if ($currentLine =~ /transcript_id (\S+)/)
    {
      $transcriptId = $1;
    }
    elsif ($currentLine =~ /transcript_id=(\S+)/)
    {
      $transcriptId = $1;
    }
    else
    {
      die "Improperly formatted gff file: transcript_id: $transcriptId";
    }
    if ($transcriptId eq $currentTranscript)
    {
      push(@currentParse, $currentLine);
    }
    else
    {
      if ($currentTranscript eq "-1")
      {
        $currentTranscript = $transcriptId;
        push(@currentParse, $currentLine);
      }
      else
      {
        my @currentParseCopy = @currentParse;
        push(@parses, \@currentParseCopy);
        @currentParse = ();
        $currentTranscript = $transcriptId;
        push(@currentParse, $currentLine);
      }
    }
  }
  push(@parses, \@currentParse);
  return \@parses;
}

#format GFF lines for AStalavista
sub format
{
  my ($inGffRef, $tempChr, $tempTranscriptId) = @_;
  my @inGff = @$inGffRef;
  my @outGff = ();
  for (my $i = 0; $i < @inGff; $i++)
  {
    my $currentLine = $inGff[$i];
    chomp($currentLine);
    if ($currentLine =~ /^#/)
    {
      next;
    }
    my @currentLineSplit = split(/\s+/, $currentLine);
    if (!($currentLineSplit[2] =~ /exon/))
    {
      die ("Found a feature that isn't an exon\n");
    }
    my ($featureStart, $featureEnd, $gffTranscriptString) = ($currentLineSplit[3], $currentLineSplit[4], $currentLineSplit[8]);
    my ($transcriptNum, $parseNum) = ($gffTranscriptString =~ /transcript_id=(\d+);parse=(\d+)/);
    my $newTranscriptLabel = "transcript_id $tempTranscriptId;";
    my $newLine = "$tempChr\t$currentLineSplit[1]\texon\t$currentLineSplit[3]\t$currentLineSplit[4]\t$currentLineSplit[5]\t$currentLineSplit[6]\t$currentLineSplit[7]\t$newTranscriptLabel\n";
    push(@outGff,$newLine);
  }
  return \@outGff;
}
