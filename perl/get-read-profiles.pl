#!/usr/bin/perl
#Get RNA-seq read counts and splice counts for a genomic region
#with the given coordinates, outputting the read and splice profiles
#as files with the given prefix.

#This script was written in part using modified code from
#three other scripts: get-gene-data.pl, make-histogram-exon.pl,
#and get_sam_reads.pl

use warnings;
use strict;


#3 parts:
#get_gene_data.pl
#make_histogram_*.pl - but just coordinates, not distances from 3' end
#get_sam_reads.pl

die "usage: get-read-profiles.pl chromosome begin end [readStrand:+-0] [transcriptStrand:+-0]  <out-prefix> <bam-file> [samtools-path]" unless @ARGV == 7 || @ARGV == 8;

my $chromosome = $ARGV[0];
my $begin = $ARGV[1];
my $end = $ARGV[2];
my $readStrand = $ARGV[3];
my $transcriptStrand = $ARGV[4];
my $outPrefix = $ARGV[5];
my $bamFile = $ARGV[6];
my $samtoolsPath = "samtools";
if (exists($ARGV[7]))
{
    $samtoolsPath = "$ARGV[7]/samtools";
}

die "readStrand must be either +, -, or 0" unless $readStrand eq "+" || $readStrand eq "-" || $readStrand eq "0"; #strand value of 0 means strand info is being ignored
die "transcriptStrand must be either +, -, or 0" unless $transcriptStrand eq "+" || $transcriptStrand eq "-" || $transcriptStrand eq "0";
die "either reads or transcript must be stranded" if $readStrand eq "0" && $transcriptStrand eq "0";

my $flag;
if ($readStrand eq "-") #see SAM format specification for an explanation of the flags
{
    $flag = "-F1023";
}
elsif ($readStrand eq "+")
{
    $flag = "-f16";
}
else
{
    $flag = "";
}

my @pileup = `$samtoolsPath mpileup -r '$chromosome:$begin-$end' $bamFile`;
my @sam = `$samtoolsPath view -h $flag $bamFile $chromosome:$begin-$end`;

my @readProfile;
my @spliceOutProfile;
my @spliceInProfile;

for (my $pCount = 0; $pCount < @pileup; $pCount++)
{
    my $currentLine = $pileup[$pCount];
    (my $coord, my $numReads, my $readString) = ($currentLine =~ /^\w+\s+(\d+)\s+\w\s+(\d+)\s+(\S+)\s+\S+/);

    my $numStrandReads = 0;
    my $numSkips = 0;
    my $length = length($readString);
    for (my $k = 0; $k < $length; $k++)
    {
        my $currentChar = substr($readString, $k, 1);
        if ($currentChar =~ /[><]/)
        {
            $numSkips++;
        }
        if ($readStrand eq "-" && $currentChar =~ /[A-Z]/ || $readStrand eq "+" && $currentChar =~ /[a-z]/ || $readStrand eq "0" && $currentChar =~ /[A-Za-z]/)
        {
            $numStrandReads++;
        }
    }

    if ($coord >= $begin && $coord <= $end)
    {
        my $internalCoord = $coord - $begin;
        $readProfile[$internalCoord] = $numStrandReads;
    }
}

for (my $sCount = 0; $sCount < @sam; $sCount++)
{
    my $currentLine = $sam[$sCount];
    if ($currentLine =~ /^@/)
    {
        next;
    }
    my @inputArray = split(/\s+/, $currentLine);
    my $readString = $inputArray[9];
    my $cigarString = $inputArray[5];
    my $startCoord = $inputArray[3];

    if ($cigarString =~ /[ID]/)
    {
       next;
    }

    my @cigarSplit = split(/[MN]/, $cigarString);

    my $cigarIndex = 0;
    my $coord = $startCoord;
    while ($cigarIndex < @cigarSplit)
    {
        my $mVal = $cigarSplit[$cigarIndex];
        $cigarIndex++;
        my $nVal = $cigarSplit[$cigarIndex];
        $cigarIndex++;
        
        if (!defined $nVal)
        {
            last;
        }
        
        $coord += $mVal; #first position with spliced reads
        if ($coord >= $begin && $coord <= $end)
        {
            my $internalCoord = $coord - $begin;
            if ($readStrand eq "-" || ($readStrand eq "0" && $transcriptStrand eq "-")) #if reverse strand, it's a splice in
            {
                if (!defined $spliceInProfile[$internalCoord]) {$spliceInProfile[$internalCoord] = 0;}
                $spliceInProfile[$internalCoord]++;
            }
            elsif ($readStrand eq "+" || ($readStrand eq "0" && $transcriptStrand eq "+")) #if forward strand, it's a splice out
            {
                if (!defined $spliceOutProfile[$internalCoord]) {$spliceOutProfile[$internalCoord] = 0;}
                $spliceOutProfile[$internalCoord]++;
            }
            else
            {
               die "Error: must have either stranded reads or stranded transcript to create splice-in and splice-out profiles";
            }
        }
        $coord += $nVal - 1; #last position with spliced reads
	if ($coord >= $begin && $coord <= $end)
	{
            my $internalCoord = $coord - $begin;
            if ($readStrand eq "-" || ($readStrand eq "0" && $transcriptStrand eq "-"))
            {
                if (!defined $spliceOutProfile[$internalCoord]) {$spliceOutProfile[$internalCoord] = 0;}
                $spliceOutProfile[$internalCoord]++;
            }
            elsif ($readStrand eq "+" || ($readStrand eq "0" && $transcriptStrand eq "+")) 
            {
                if (!defined $spliceInProfile[$internalCoord]) {$spliceInProfile[$internalCoord] = 0;}
                $spliceInProfile[$internalCoord]++;
            }
            else
            {
               die "Error: must have either stranded reads or stranded transcript to create splice-in and splice-out profiles";
            }
        }
    }
}

open(READ_PROFILE, "> $outPrefix.readProfile");
open(SPLICE_IN_PROFILE, "> $outPrefix.spliceInProfile");
open(SPLICE_OUT_PROFILE, "> $outPrefix.spliceOutProfile");

print(READ_PROFILE "@ $chromosome $begin $end\n");
print(SPLICE_OUT_PROFILE "@ $chromosome $begin $end\n");
print(SPLICE_IN_PROFILE "@ $chromosome $begin $end\n");

for (my $m = 0; $m < $end - $begin + 1; $m++)
{
    if (!defined $readProfile[$m])
    {
        $readProfile[$m] = 0;
    }
    print(READ_PROFILE "$readProfile[$m]\n");
}

for (my $n = 0; $n < $end - $begin + 1; $n++)
{
    if (!defined $spliceOutProfile[$n])
    {
        $spliceOutProfile[$n] = 0;
    }
    print(SPLICE_OUT_PROFILE "$spliceOutProfile[$n]\n");
}

for (my $o = 0; $o < $end - $begin + 1; $o++)
{
    if (!defined $spliceInProfile[$o])
    {
        $spliceInProfile[$o] = 0;
    }
    print(SPLICE_IN_PROFILE "$spliceInProfile[$o]\n");
}

