#!/usr/bin/perl
#Compare predictions to GFF annotations using AStalavista and check support
#
#See the comment starting at line 98 for information about repurposing this
#script to perform other tasks using pairwise comparison

use warnings;
use strict;
use Getopt::Std;

sub format;
sub compare;
sub splitIntoParses;
sub getNucleotides;

our $opt_s;
getopts('s:');

die "compare-predictions-asta-support.pl [-s samtoolspath] <annotation-gff> <prediction-gff> <bamfile> <astalavista-path>\n" unless @ARGV==4;

my $samtoolsPath = "samtools";
if (defined($opt_s))
{
  $samtoolsPath = "$opt_s/samtools";
}

my ($aFile, $pFile, $bamfile, $astapath) = @ARGV;

open(A_FILE, $aFile);
open(P_FILE, $pFile);

my @annotation = (<A_FILE>);
my @prediction = (<P_FILE>);

#Split the annotation and prediction into arrays of parses
#Each parse is itself an array of GFF lines
my $aSplitRef = splitIntoParses(\@annotation);
my $pSplitRef = splitIntoParses(\@prediction);

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
  my $endBound = ($pFirstLineEnd > $pLastLineEnd ? $pFirstLineEnd : $pLastLineEnd);
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
    my $pBegin = ($pFirstLineBegin < $pLastLineBegin ? $pFirstLineBegin : $pLastLineBegin);
    my $aBegin = ($aFirstLineBegin < $aLastLineBegin ? $aFirstLineBegin : $aLastLineBegin);
    my $offset = ($pBegin < $aBegin ? $pBegin : $aBegin);
    #get the gene structure aways and find the number of differences
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

  #Here you could perform any other task with the current prediction and
  #its closest annotation (such as comparing proteins). To do so, use the
  #$currentPred and $minDiffAnn variables, references to arrays containing
  #the GFF lines of the prediction and its closest annotation

  #format the prediction and annotation to prepare them for comparison
  my $afRef = &format($minDiffAnn, "1");
  my $pfRef = &format($currentPred, "2");
  print("$currentPredID\t$currentAnnID\t");
  #compare the prediction with the closest annotation using AStalavista
  compareAsta($afRef, $pfRef);
  #compareEnds($afRef, $pfRef);
}

#Use AStalavista to find the alternative splice events in two GFF files
#and compare their support
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
    print("no difference\n");
  }
  elsif (@structureList > 1)
  {
    print("more than one AS event\n");
  }
  else
  {
    my @structure = @{$structureList[0]};
    checkSupport(\@structure);
  }
}

#Figure out what type of AS event is occurring, then use samtools to see
#how many reads support the annotated and predicted versions of the gene
sub checkSupport
{
  my $a = 8; #required number of bases the spliced read must have on each side of the splice site
  my ($struct) = @_;
  my @structure = @$struct;
  my ($start, $end, $scStart, $scEnd, $code, $chr, $trA, $trB, $strand) = @structure;
  #print("$start $end $scStart $scEnd $code $chr\n");
  print("$code\t");
  my ($samStart, $samEnd);
  if ($code eq "1^,2^") #alt donor
  {
    $samStart = ($scStart < $scEnd ? $scStart : $scEnd) - 92;
    $samEnd = $end + 92;
  }
  elsif ($code eq "1-,2-") #alt acceptor
  {
    $samStart = $start - 92;
    $samEnd = ($scStart > $scEnd? $scStart : $scEnd) + 92;
  }
  elsif ($code eq "1-2^,0" || $code eq "0,1-2^") #exon skipping
  {
    $samStart = $start;
    $samEnd = $end;
  }
  else
  {
    print("unable to check support: unsupported structure\n");
    return;
  }
  my @samreads = `$samtoolsPath view $bamfile $chr:$samStart-$samEnd`;
  my ($countA, $countB) = (0,0);
  for (my $i = 0; $i < @samreads; $i++)
  {
    my $read = $samreads[$i];
    my @readSplit = split(/\s+/, $read);
    my ($rnaStart, $cigar) = ($readSplit[3], $readSplit[5]);
    #length calculations: generally all coordinates are beginnings/ends of exon,
    #so if calculating the length of an intron between an end and a start, subtract
    #one; if calculating the length of an exon between a start and an end, add one
    if (($code eq "1^,2^" && $strand eq "+") || ($code eq "1-,2-" && $strand eq "-")) #forward alt donor, reverse alt acceptor
    {
      my $junctionLengthA = $end - $scStart - 1;
      my $junctionLengthB = $end - $scEnd - 1;
      #check for position (>=8bp on each side)
      if ($cigar =~ /(\d+)M${junctionLengthA}N(\d+)M/ && $1 >= $a && $2 >= $a && ($rnaStart + $1 - 1 == $scStart))
      {
        $countA++;
      }
      elsif ($cigar =~ /(\d+)M${junctionLengthB}N(\d+)M/ && $1 >= $a && $2 >= $a && ($rnaStart + $1 - 1 == $scEnd))
      {
        $countB++;
      }
    }
    elsif (($code eq "1-,2-" && $strand eq "+") || ($code eq "1^,2^" && $strand eq "-")) #forward alt acceptor, reverse alt donor
    {
      my $junctionLengthA = $scStart - $start - 1;
      my $junctionLengthB = $scEnd - $start - 1;
      if ($cigar =~ /(\d+)M${junctionLengthA}N(\d+)M/ && $1 >= $a && $2 >= $a && ($rnaStart + $1 - 1 == $start)) #same as above
      {
        $countA++;
      }
      elsif ($cigar =~ /(\d+)M${junctionLengthB}N(\d+)M/ && $1 >= $a && $2 >= $a && ($rnaStart + $1 - 1 == $start))
      {
        $countB++;
      }
    }
    elsif ($code eq "1-2^,0" || $code eq "0,1-2^") #exon skipping
    {
      my ($internalExonLength, $littleIntronALength, $littleIntronBLength, $bigIntronLength);
      if ($strand eq "+")
      {
        $internalExonLength = $scEnd - $scStart + 1;
        $littleIntronALength = $scStart - $start - 1;
        $littleIntronBLength = $end - $scEnd - 1;
        $bigIntronLength = $end - $start - 1;
      }
      elsif ($strand eq "-")
      {
        $internalExonLength = $scStart - $scEnd + 1;
        $littleIntronALength = $scEnd - $start - 1;
        $littleIntronBLength = $end - $scStart - 1;
        $bigIntronLength = $end - $start - 1;
      }
      if ($cigar =~ /(\d+)M${littleIntronALength}N${internalExonLength}M${littleIntronBLength}N(\d+)M/
        && ($rnaStart + $1 - 1 == $start))
      {
        if ($1 >= $a && $2 >= $a) #check here so that the CIGAR string doesn't rematch a later regex
        {
          $countB++;
        }
      }
      elsif ($cigar =~ /(\d+)M${littleIntronALength}N(\d+)M/ && $1 >= $a && $2 >= $a
        && $2 <= $internalExonLength)
      {
        if ($rnaStart + $1 - 1 != $start) #looks like x1Mx2Nx3Mx4Nx5M, where the end of the read overlaps with little intron A
        {
          if ($cigar =~ /(\d+)M(\d+)N(\d+)M(\d+)N(\d+)M/)
          {
            my ($l1, $l2, $l3, $l4, $l5) = ($1, $2, $3, $4, $5);
            if ($rnaStart + $l1 + $l2 + $l3 - 1 == $start)
            {
              $countB++;
            }
          }
        }
        else
        {
          $countB++;
        }
      }
      elsif ($cigar =~ /(\d+)M${littleIntronBLength}N(\d+)M/ && $1 >= $a && $2 >= $a
        && $1 <= $internalExonLength && ($rnaStart + $1 - 1 == ($strand eq "+" ? $scEnd : $scStart)))
      {
        $countB++;
      }
      elsif ($cigar =~ /(\d+)M:${bigIntronLength}N:(\d+)M/ && $1 >= $a && $2 >= $a && ($rnaStart + $1 - 1 == $start))
      {
        $countA++;
      }
    }
  }
  my ($countAnn, $countPred);
  if ($trA == 1)
  {
    ($countAnn, $countPred) = ($countA, $countB);
  }
  else
  {
    ($countAnn, $countPred) = ($countB, $countA);
  }
  print("$countPred\t$countAnn\n");
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
  print("Ends: ");
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
  print("\n");
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
    my @linesplit = split(/\s+/, $astafile[$i]);
    my $astaCode;
    my ($trA, $trB);
    my ($scStart, $scEnd);
    my ($start, $end) = ($linesplit[3], $linesplit[4]);
    my $chr = $linesplit[0];
    my $strand = $linesplit[6];
    if ($astafile[$i] =~ /structure "(\S+)"/)
    {
      $astaCode = $1;
    }
    else
    {
      print("Error in AStalavista output pattern match\n");
    }
    if ($astafile[$i] =~ /splice_chain "[,\-\^]*(\d+)[,\-\^]*(\d+)[,\-\^]*\S*"/)
    {
      ($scStart, $scEnd) = ($1, $2);
    }
    else
    {
      print("Error in AStalavista output pattern match\n");
      print("$astafile[$i]\n");
    }
    if ($astafile[$i] =~ /transcript_id "(\d+);,(\d+);"/)
    {
      ($trA, $trB) = ($1, $2);
    }
    else
    {
      print("Error in AStalavista output pattern match\n");
    }
    my @structure = ($start, $end, $scStart, $scEnd, $astaCode, $chr, $trA, $trB, $strand);
    push(@structureList, \@structure);
  }
  return \@structureList;
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

#Format GFF lines for AStalavista
sub format
{
  my ($inGffRef, $tempTranscriptId) = @_;
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
    my ($featureStart, $featureEnd) = ($currentLineSplit[3], $currentLineSplit[4]);
    my $newTranscriptLabel = "transcript_id $tempTranscriptId;";
    my $newLine = "$currentLineSplit[0]\t$currentLineSplit[1]\texon\t$currentLineSplit[3]\t$currentLineSplit[4]\t$currentLineSplit[5]\t$currentLineSplit[6]\t$currentLineSplit[7]\t$newTranscriptLabel\n";
    push(@outGff,$newLine);
  }
  return \@outGff;
}
