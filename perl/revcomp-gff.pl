#!/usr/bin/perl
#Reverse complement a gff prediction output by Decoder

use strict;
use warnings;

die "revcomp-gff.pl  <gff-file>" unless @ARGV == 1;

my $gffFile = "$ARGV[0]";

open(GFF, $gffFile);
my @newgff;
my @gff = <GFF>;
close (GFF);

my @topComments = ();
my @parses = ();
my @sequenceLengthComment = ();
my @scores = ();


#get transcript length
my $transcriptLength = -1;

for (my $i = 0; $i < @gff; $i++)
{
    if ($gff[$i] =~ /# sequence length: (\d+)/)
    {
        $transcriptLength = $1;
    }
}

#get each parse, put it in an array, and put the ref in @parses
my $currentLineIndex = 0;

while ($currentLineIndex < @gff)
{
    my $currentLine = $gff[$currentLineIndex];
    chomp($currentLine);
    if ($currentLine =~ /##/)
    {
        push(@topComments, $currentLine);
        $currentLineIndex++;
        next;
    }
    if ($currentLine =~ /# sequence length:/)
    {
        push(@sequenceLengthComment, $currentLine);
        $currentLineIndex++;
        next;
    }
    if ($currentLine =~ /score/)
    {
        push(@scores, $currentLine);
        $currentLineIndex++;
        next;
    }
    my ($parseNum) = ($currentLine =~ /parse=(\d+);/);
    my $parseIndex = (defined($parseNum) ? $parseNum - 1 : 0);
    
    if (!defined($parses[$parseIndex]))
    {
        my @newArray = ();
        $parses[$parseIndex] = \@newArray;
    }

    push(@{$parses[$parseIndex]}, $currentLine);

    $currentLineIndex++;
}


my @revcompParses = ();
for (my $i = 0; $i < @parses; $i++)
{
    $revcompParses[$i] = revcompParse($parses[$i], $transcriptLength);
}

printElements(\@topComments);
for (my $i = 0; $i < @revcompParses; $i++)
{
    printElements($revcompParses[$i]);
}

printElements(\@sequenceLengthComment);
printElements(\@scores);

##########################
#PRINT AN ARRAY'S ELEMENTS,
#EACH ON ONE LINE
##########################
sub printElements
{
    my ($ref) = @_;
    my @list = @{$ref};
    for (my $i = 0; $i < @list; $i++)
    {
        my $currentLine = $list[$i];
        chomp($currentLine);
        print("$currentLine\n");
    }
}

#####################
#REVERSE AN ARBITRARY
#LIST
#####################
sub revList
{
    my ($listRef) = @_;
    my @list = @{$listRef};
    my @newList = ();
    for (my $i = 0; $i < @list; $i++)
    {
        $newList[(@list - 1) - $i] = $list[$i];
    }
    return \@newList;
}

######################
#REVERSE COMPLEMENT AN
#INDIVIDUAL PARSE
######################

sub revcompParse
{
my ($parseRef, $transcriptLength) = @_;
my @parse = @{$parseRef};
my @newParse = ();
my $length = $transcriptLength; #unnecessary duplicate?
my ($numFeatures, $maxTranscriptNum) = (0,0);
for (my $h = 0; $h < @parse; $h++) #do an initial loop to figure out max transcript num, num features
{
    my $currentLine = $parse[$h];
    chomp($currentLine);
    my @lineSplit = split(/\s+/, $currentLine);
    my $TPString = $lineSplit[8];
    my @TPSplit = split(/;/, $lineSplit[8]);
    my $transcriptString = $TPSplit[0];
    my ($transcriptNum) = ($transcriptString =~ /(\d+)/);
    if ($transcriptNum > $maxTranscriptNum)
    {
        $maxTranscriptNum = $transcriptNum;
    }
    $numFeatures++;
}

for (my $i = 0; $i < @parse; $i++)
{
    my $currentLine = $parse[$i];
    chomp($currentLine);

    my @lineSplit = split(/\s+/, $currentLine);
    my ($featureStart, $featureEnd) = ($lineSplit[3], $lineSplit[4]);
    my ($intFeatureStart, $intFeatureEnd) = ($featureStart - 1, $featureEnd - 1); #convert to 0-based
    my $featureStrand = $lineSplit[6];
    my $oldTPString = $lineSplit[8];
    my @oldTPSplit = split(/;/, $oldTPString);
    my $oldTranscriptString = "$oldTPSplit[0];";

    my $oldParseString = (defined($oldTPSplit[1]) ? "$oldTPSplit[1];" : "");
    my ($newIntFeatureStart, $newIntFeatureEnd) = ($length - 1 - $intFeatureEnd, $length - 1 - $intFeatureStart);
    my ($newFeatureStart, $newFeatureEnd) = ($newIntFeatureStart + 1, $newIntFeatureEnd + 1); #convert back to 1-based
    my $newFeatureStrand = ($featureStrand eq "+") ? "-" : "+";    

    my ($oldTranscriptNum) = ($oldTranscriptString =~ /(\d+)/);
    my $newTranscriptNum = $maxTranscriptNum + 1 - $oldTranscriptNum;
  
    my $newTranscriptString = "transcript_id=${newTranscriptNum};";
    my $newParseString = $oldParseString;

    my $newLine = "$lineSplit[0]\t$lineSplit[1]\t$lineSplit[2]\t$newFeatureStart\t$newFeatureEnd\t$lineSplit[5]\t$newFeatureStrand\t$lineSplit[7]\t$newTranscriptString$newParseString";

    if (defined($lineSplit[9]))
    {
      $newLine = "$newLine $lineSplit[9]";
    }

    $newParse[(@parse - 1) - $i] = $newLine; #reverse the elements here
}
my $newParseRef = \@newParse;
return $newParseRef;
}

#####
#WhT?
#####
