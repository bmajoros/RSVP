#!/usr/bin/perl
#Compare predictions to GFF annotations and produce pdiff output
#
#See the comment starting at line 98 for information about repurposing this
#script to perform other tasks using pairwise comparison
#
#Comment out the system call on line 135 if you want to save the pdiff files

use warnings;
use strict;

sub format;
sub splitIntoParses;
sub getNucleotides;


die "compare-predictions-pdiff.pl <annotation-gff> <prediction-gff> <pdiff-out-dir> <script-dir>\n" unless @ARGV==4;

my ($aFile, $pFile, $outDir, $scriptDir) = @ARGV;

open(A_FILE, $aFile);
open(P_FILE, $pFile);

my ($chunk) = ($aFile =~ /(\w+).gff/) or die "annotation GFF file name error: expecting <chunk-name>.gff";

my @annotation = (<A_FILE>);
my @prediction = (<P_FILE>);

#Split the annotation and prediction into arrays of parses
#Each parse is itself an array of GFF lines
my $aSplitRef = splitIntoParses(\@annotation);
my $pSplitRef = splitIntoParses(\@prediction);

my $prefix; #used to concatenate multiple pdiff files into one file of SAM-type strings

#Compare each prediction to all of the annotations and find the closest match
for (my $i = 0; $i < @$pSplitRef; $i++)
{
  my $currentPred = $pSplitRef->[$i];
  #find the end of the prediction, since we only compare up to that point
  my $pFirstLine = $currentPred->[0];
  if (!(defined($pFirstLine))) {die "error: prediction has no lines. Corresponding annotation file: $aFile"}
  my @pFirstLineSplit = split(/\s+/,$pFirstLine);
  my ($pFirstLineBegin, $pFirstLineEnd) = ($pFirstLineSplit[3], $pFirstLineSplit[4]);
  my $pLastLine = $currentPred->[@$currentPred-1];
  my @pLastLineSplit = split(/\s+/,$pLastLine);
  my ($pLastLineBegin, $pLastLineEnd) = ($pLastLineSplit[3], $pLastLineSplit[4]);
  my $pBegin = ($pFirstLineBegin < $pLastLineBegin ? $pFirstLineBegin : $pLastLineBegin);
  my $pEnd = ($pFirstLineEnd > $pLastLineEnd ? $pFirstLineEnd : $pLastLineEnd);
  my $predStrand = $pFirstLineSplit[6];
  #compare pairwise to all the annotations
  my $minDiff = 9**9**9; #infinity
  my $minDiffAnn;
  my $minDiffAnnNum;
  for (my $j = 0; $j < @$aSplitRef; $j++)
  {
    my $currentAnn = $aSplitRef->[$j];
    #find starting coordinate for gene structure arrays
    my $aFirstLine = $currentAnn->[0];
    my @aFirstLineSplit = split(/\s+/,$aFirstLine);
    my ($aFirstLineBegin, $aFirstLineEnd) = ($aFirstLineSplit[3], $aFirstLineSplit[4]);
    my $aLastLine = $currentAnn->[@$currentAnn-1];
    my @aLastLineSplit = split(/\s+/,$aLastLine);
    my ($aLastLineBegin, $aLastLineEnd) = ($aLastLineSplit[3], $aLastLineSplit[4]);
    my $aBegin = ($aFirstLineBegin < $aLastLineBegin ? $aFirstLineBegin : $aLastLineBegin);
    my $offset = ($pBegin < $aBegin ? $pBegin : $aBegin);
    my $annStrand = $aFirstLineSplit[6];
    die "error: annotation and prediction strands do not match" if !($annStrand eq $predStrand);
    #get the gene structure arrays and find the number of differences
    my $predNucleotides = getNucleotides($currentPred, $offset);
    my $annNucleotides = getNucleotides($currentAnn, $offset);
    my $differences = 0;
    for (my $k = 0; $k <= $pEnd - $offset; $k++)
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
  my ($currentPredID, $currentAnnID);
  if ($pFirstLine =~ /transcript_id[\s=]\"*(\w*\.*\d+\.*\d*)\"*;/)
  {
    $currentPredID = $1;
  }
  else
  {
    die("Prediction transcript ID formatted incorrectly");
  }
  if ($minDiffAnn->[0] =~ /transcript_id[\s=]\"*(\w*\.*\d+)\"*;/)
  {
    $currentAnnID = $1;
  }
  else
  {
    die("Annotation transcript ID formatted incorrectly");
  }

  #Here you could perform any other task with the current prediction and
  #its closest annotation (such as comparing proteins). To do so, use the
  #$currentPred and $minDiffAnn variables, references to arrays containing
  #the GFF lines of the prediction and its closest annotation
  my $temp1 = "temp1.gff";
  my $temp2 = "temp2.gff";
  my ($predTranscriptNum, $parseNum) = ($currentPredID =~ /\w+\.(\d+)\.(\d+)/) or die "error in prediction transcript ID formatting";
  my ($gene, $annTranscriptNum) = ($currentAnnID =~ /(\w+)\.(\d+)/) or die "error in annotation transcript ID formatting";
  my $out = "$outDir/$gene.${annTranscriptNum}_$predTranscriptNum.$parseNum.pdiff";
  printGFF($minDiffAnn, $temp1);
  printGFF($currentPred, $temp2);
  makePdiff($temp1, $temp2, $out);
  my $header = "#annotation: $currentAnnID\tprediction: $currentPredID\tstrand: $predStrand\tchunk: $chunk";
  makeHeader($out, $header);
  system("rm $temp1");
  system("rm $temp2");
  if (!(defined($prefix)))
  {  
    $prefix = "$outDir/$gene";
  }
  else
  {
    die "all annotations do not have the same gene name" if !($prefix eq "$outDir/$gene");
  }
}
makeStrings($prefix);
system("rm $prefix*.pdiff");

sub printGFF
{
  my ($ref, $file) = @_;
  open(GFF, "> $file");
  for (my $i = 0; $i < @$ref; $i++)
  {
    print(GFF $ref->[$i]);
  }
}

sub makePdiff
{
  my ($ann, $pred, $out) = @_;
  system("perl $scriptDir/make-pdiff.pl $ann $pred > $out");
}

sub makeStrings
{
  my $prefix = shift;
  system("perl $scriptDir/process-pdiffs.pl $prefix $scriptDir");
}

sub makeHeader
{
  my ($out, $header) = @_;
  open(OLD_PDIFF, $out);
  my @oldPdiff = <OLD_PDIFF>;
  close(OLD_PDIFF);
  open(NEW_PDIFF, "> $out");
  print(NEW_PDIFF "$header\n");
  print(NEW_PDIFF @oldPdiff);
  close(NEW_PDIFF);
}

#Make an array for each parse/annotation where array[i] is 1 if i is in
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

#Return an array of parses, each itself an array of lines
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
    if ($currentLine =~ /transcript_id \"*(\S+)\"*/)
    {
      $transcriptId = $1;
    }
    elsif ($currentLine =~ /transcript_id=\"*(\S+)\"*/)
    {
      $transcriptId = $1;
    }
    else
    {
      die "Improperly formatted gff file: $currentLine";
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
