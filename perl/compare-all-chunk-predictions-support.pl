#!/usr/bin/perl
#Run compare-predictions-asta-support.pl on all chunks in a folder

use warnings;
use strict;
use Getopt::Std;

our $opt_s;
getopts('s:');

die "compare-all-chunk-predictions-support.pl [-s samtoolspath] <annotation-dir> <prediction-dir> <bam-file> <astalavista-path> <script-dir>" unless @ARGV==5;

print("pred.#\t\tann.#\t\tdiff\tRS new\tRS old\n");

my ($dataDir, $outDir, $bamFile, $astaPath, $scriptDir) = @ARGV;

my $files = `ls $dataDir`;
my @files = split(/\s+/, $files);

for my $file (@files)
{
  if ($file =~ /(\d+).gff/)
  {
    my $fileA = "$dataDir/$1.gff";
    my $fileB = "$outDir/$1.gff";
    system("perl $scriptDir/keep-center-transcripts.pl $fileB > cacp-temp.gff");
    if (defined($opt_s))
    {
      system("perl $scriptDir/compare-predictions-asta-support.pl -s $opt_s $fileA cacp-temp.gff $bamFile $astaPath");
    }
    else
    {
      system("perl $scriptDir/compare-predictions-asta-support.pl $fileA cacp-temp.gff $bamFile $astaPath");
    }
    system("rm cacp-temp.gff");
  }
}
