#!/usr/bin/perl
#Calculates the exon score histogram (array of Bernoulli distributions)
#for each distance d away from the stop codon. Divides the number of
#genes with reads >= threshold at a given distance by the number of
#genes extending at least that far. Throws away genes with evidence of
#alternative splicing (reference skips inside exon or UTR or stranded 
#reads inside an intron). The first line of the output is a comment with
#the sample size (i.e. number of genes not thrown away); the rest of the
#output is distance-probability pairs.

use strict;
use warnings;

die "usage: make-histogram-exon.pl threshold <data-directory> [-s]" unless @ARGV == 2 || @ARGV == 3;

my $threshold = $ARGV[0]; #check for at least this many reads
my $dataDir = $ARGV[1];
my $stranded = (defined($ARGV[2]) && $ARGV[2] eq "-s" ? 1 : 0);

my $chunkFile = "single-transcript-genes.gff";

open(CHUNK_FILE, "$dataDir/$chunkFile") || die "failed to open chunk file\n";
my @chunkInfo = <CHUNK_FILE>;
my $i = 0;

my $numFiles = @chunkInfo;
my $numSkipped = 0;

my @globalHistogram;
my @lengthHistogram;

GENE: while ($i < $numFiles)
{
    (my $chunkStart, my $chunkEnd, my $chunkID) = ($chunkInfo[$i] =~ /^\w+\s+\w+\s+\w+-*\w*\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+[\+-]\s+[\w\.]\s+\w+=[\w\.]+;?\s+\w+=(\d+);?/);
    
    open(PILEUP_FILE, "$dataDir/$chunkID.pileup") || print "failed to open pileup file $chunkID\n";
    
    open(GFF_FILE, "$dataDir/$chunkID.gff") || print "failed to open GFF file $chunkID\n";

    #read an annotation to get a list of exon information
    my @exonStartList;
    my @exonEndList;
    my @exonDistanceList;
    my $exonCounter = 0;
    my $geneStrand;

    while (<GFF_FILE>)
    {
	my $currentLine = $_;
	(my $exonStart, my $exonEnd, my $strand) = ($currentLine =~ /^-?\w+\s+\w+\s+\w+-\w+\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+([\+-])\s+\w\s+\w+\s+[""\w\.]+/);
	if (!defined $geneStrand)
	{
	    $geneStrand = $strand;
	}
	$exonStartList[$exonCounter] = $exonStart;
        $exonEndList[$exonCounter] = $exonEnd;
        $exonCounter++;
    }
    #end annotation gff read
    
    my $numExons = @exonStartList;

    #calculate distances from 5' end of exon to 3' end of strand
    for (my $distCounter = $numExons - 1; $distCounter >= 0; $distCounter--)
    {
	my $curExonLength = $exonEndList[$distCounter] - $exonStartList[$distCounter] + 1;
	if ($distCounter >= $numExons - 1)
	{
	    $exonDistanceList[$distCounter] = $curExonLength - 1;
	}
	else
	{
	    $exonDistanceList[$distCounter] = $curExonLength + $exonDistanceList[$distCounter + 1]; #+1-1 = 0
	}
    }
    
    #build histogram using pileup data
    my @currentHistogram;
    LINE: while (<PILEUP_FILE>)
    {
	(my $coord, my $numReads, my $readString) = /^\w+\s+(\d+)\s+\w\s+(\d+)\s+(\S+)\s+\S+/;

	my $numStrandReads = 0;
	my $numStrandSkips = 0;
        my $length = length($readString);
        for (my $k = 0; $k < $length; $k++)
        {
            my $currentChar = substr($readString, $k, 1);
            if ($geneStrand eq "+" && $currentChar eq "<" || $geneStrand eq "-" && $currentChar eq ">" || $stranded == 0 && $currentChar =~ /[><]/)
            {
                $numStrandSkips++;
            }
            if ($geneStrand eq "+" && $currentChar =~ /[a-z]/ || $geneStrand eq "-" && $currentChar =~ /[A-Z]/ || $stranded == 0 && $currentChar =~ /[A-Za-z]/)
            {
                $numStrandReads++;
            }
        }
	
	my $normalizedCoord = $coord; #all-genes exon gffs don't have normalized coordinates
	my $distIntoExon = -1;
	
	my $currentExonDistance = -1;
	my $currentExon = -1;
	my $currentState = ""; #exon, intron, or UTR
	my $foundPos = 0;
	
	
	my $j = 0;
	my $increment = 1;

	if ($geneStrand eq "+")
	{
	    $j = 0;
	    $increment = 1;
	}
	else
	{
	    $j = $numExons - 1;
	    $increment = -1;
	}


	if ($normalizedCoord < $exonStartList[$j])
	{
	    $currentState = "UTR";
	    $foundPos = 1;
	}
	if ($normalizedCoord >= $exonStartList[$j] && $normalizedCoord <= $exonEndList[$j])
	{
	    $currentState = "exon";
	    $currentExon = $j;
	    $foundPos = 1;
	}
	
	while ($j < $numExons && $j > -1 && $foundPos == 0)
	{
	    if ($normalizedCoord < $exonStartList[$j] && $normalizedCoord > $exonEndList[$j - 1])
	    {
		$currentState = "intron";
		$foundPos = 1;
	    }
	    if ($normalizedCoord >= $exonStartList[$j] && $normalizedCoord <= $exonEndList[$j])
	    {
		$currentState = "exon";
		$currentExon = $j;
		$foundPos = 1;
	    }
	    $j += $increment;
	}
	if ($foundPos == 0)
	{
	    $currentState = "UTR";
	    $foundPos = 1;
	}

        if ($currentState eq "UTR")
        {
            if ($numStrandSkips > 0)
            {
                $i++;
                $numSkipped++;
                next GENE;
            }            
            else
            {
                next LINE;
            }
        }
        if ($currentState eq "intron")
        {
            if ($numStrandReads > 0)
            {
                $i++;
                $numSkipped++;
                next GENE;
            }
            else
            {
                next LINE;
            }
        }
        if ($currentState eq "exon")
        {    
            if ($numStrandSkips > 0)
            {
                $i++;
                $numSkipped++;
                next GENE;
            }
        }

	$currentExonDistance = $exonDistanceList[$currentExon];
	if ($geneStrand eq "+")
	{
	    $distIntoExon = $normalizedCoord - $exonStartList[$currentExon];
	}
	else
	{
	    $distIntoExon = $exonEndList[$currentExon] - $normalizedCoord;
	}

	my $dist3End = $currentExonDistance - $distIntoExon;
	$currentHistogram[$dist3End] = $numStrandReads;
    }
    
    #add to the big histogram
    my $curHistSize = @currentHistogram;
    for (my $b = 0; $b < $curHistSize; $b++)
    {
	if (!defined $currentHistogram[$b])
	{
	    $currentHistogram[$b] = 0;
	}
	
	if (!defined $globalHistogram[$b])
	{
	    $globalHistogram[$b] = 0;
	}
	
	if ($currentHistogram[$b] >= $threshold)
	{
	    $globalHistogram[$b]++;
	}
    }

    #calculate length
    my $currentLength = 0;
    for (my $e = 0; $e < @exonStartList; $e++)
    {
        $currentLength += ($exonEndList[$e] - $exonStartList[$e] + 1);
    }
    $currentLength -= 1; #we want the length as the maximum distance from 3' end

    #add to the big length histogram
    for (my $f = 0; $f <= $currentLength; $f++)
    {
        if (!defined $lengthHistogram[$f])
        {
            $lengthHistogram[$f] = 0;
        }
        $lengthHistogram[$f]++;
    }

    $i++;
}

$numFiles -= $numSkipped;
print("#number of genes used: $numFiles\n");

#print big histogram
my $globalHistSize = @globalHistogram;

for (my $c = 0; $c < $globalHistSize; $c++)
{
    my $probability = $globalHistogram[$c] / $lengthHistogram[$c];
    print("$c $probability\n");
}



#chunk-gff-regex: ^\w+\s+\w+\s+\w+-*\w*\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+[\+-]\s+[\w\.]\s+\w+=[\w\.]+;?\s+\w+=(\d+);?
#OLD chunk-gff-regex: ^\w+\s+\w+\s+\w+-*\w*\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+[\+-]\s+\w\s+\w+=[\w\.]+;\s+\w+=(\d+)
#pileup-regex: ^\w+\s+\d+\s+\w\s+(\d+)\s+\S+\s+\S+
#OLD other-gff-regex: ^-?\w+\s+\w+\s+\w+-\w+\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+([\+-])\s+\w\s+\w+\s+["\w\.]+
#other-gff-regex: ^-?\w+\s+\w+\s+\w+-\w+\s+(\d+)\s+(\d+)\s+-?\w*\.?\w*\s+([\+-])\s+\w\s+\w+\s+[""\w\.]+
