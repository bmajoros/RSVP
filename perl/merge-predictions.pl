#!/usr/bin/perl
#Merge forward and reverse strand predictions. Right now this script simply
#picks the strand with more rna-seq reads

use warnings;
use strict;

die "merge-predictions.pl <out-prefix> <data-prefix>" unless @ARGV == 2;

my ($outPrefix, $dataPrefix) = @ARGV;

my ($reverseProfile, $forwardProfile) = ("$dataPrefix.readProfile", "$dataPrefix.forward.readProfile");
my ($reverseCount, $forwardCount) = (countReads($reverseProfile), countReads($forwardProfile));

my @result;
my ($reverseGFF, $forwardGFF) = ("$outPrefix.reverse.gff", "$outPrefix.forward.gff");
open(REVERSE_GFF, $reverseGFF);
open(FORWARD_GFF, $forwardGFF);

my @reverse = <REVERSE_GFF>;
my @forward = <FORWARD_GFF>;
if ($reverseCount > $forwardCount)
{
    @result = @reverse;
}
else
{
    @result = @forward;
}

#TEMP
#my ($chunk) = ($reversePath =~ /(\d+).reverse.gff$/);
#open(COORDS, "/home/ohler/nl53/rnaseq/niel_decoder_runs/chunks2_test/$chunk.coords.gff");
#my $coords = <COORDS>;
#chomp($coords);
#my @coordsSplit = split(/\s+/, $coords);
#my $strand = $coordsSplit[6];
#if ($strand eq "-")
#{
#    @result = @reverse;
#}
#else
#{
#    @result = @forward;
#}
#end TEMP

for (my $k = 0;  $k < @result; $k++)
{
    print($result[$k]);
}

sub countReads
{
    my $fileName = shift(@_);
    open(PROFILE, $fileName) or die "unable to read profile $fileName";
    my @profile = <PROFILE>;
    my $count = 0;
    for (my $i = 1; $i < @profile; $i++)
    {
        my $currentNum = $profile[$i];
        chomp($currentNum);
        $count += $currentNum;
    }
    return $count;
}
