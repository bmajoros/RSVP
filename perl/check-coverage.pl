#!/usr/bin/perl
#A different way of eliminating genes with no RNA-seq reads that fits with
#the new pipeline. Removes genes with less than the given number of
#reads per kilobase. 
use warnings;
use strict;
use Getopt::Std;

our $opt_s;
getopts('s:');

die "usage: check-coverage.pl [-s samtoolspath] <chunk-dir> <threshold> <bam-file>" unless @ARGV == 3;

my $samtoolsPath = "samtools";
if (defined($opt_s))
{
  $samtoolsPath = "$opt_s/samtools";
}
my $bamFile = $ARGV[2];

my $dir = $ARGV[0];
my $files = `ls $dir/*.coords.gff`;
my @files = split(/\s+/, $files);
my $numZero = 0;
my $threshold = $ARGV[1];
my $numSubThreshold = 0;

for my $file (@files)
{
    open(COORDS, $file) || die "couldn't open file $file";
    my $flag = 0;
    (my $fileNum) = ($file =~ /$dir\/(\d+).coords.gff/);
    my $coords = <COORDS>;
    chomp($coords);
    my ($chr, $set, $feature, $start, $end, $score, $strand, $placeholder) = split(/\s+/, $coords);
    #set, feature, score, and placeholder are all currently not used

    my $count = `$samtoolsPath view $bamFile '$chr:$start-$end' | wc`;
    my @countSplit = split(/\s+/, $count);
    my $numReads = $countSplit[1];
    my $length = $end - $start + 1;
    my $rpkb = $numReads / ($length / 1000);
    if ($rpkb == 0) {$numZero++;}
    if ($rpkb < $threshold) {$numSubThreshold++;}

    if($rpkb < $threshold)
    {
        system("rm $dir/$fileNum.*");
    }
}
#print("num files with no reads: $numZero\n");
#print("num files with rpkb < $threshold: $numSubThreshold\n");
