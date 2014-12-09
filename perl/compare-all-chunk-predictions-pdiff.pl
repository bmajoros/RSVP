#!/usr/bin/perl
use strict;

#Run compare-predictions-pdiff.pl on all chunks in a folder

die "compare-all-chunk-predictions-pdiff.pl <annotation-dir> <prediction-dir> <pdiff-out-dir> <script-dir>" unless @ARGV==4;

my ($annDir, $predDir, $outDir, $scriptDir) = @ARGV;

my $files = `ls $annDir`;
my @files = split(/\s+/, $files);

for my $file (@files)
{
  if ($file =~ /(\d+).gff/)
  {
    my $fileA = "$annDir/$1.gff";
    my $fileB = "$predDir/$1.gff";
    system("perl $scriptDir/keep-center-transcripts.pl $fileB > cacpp-temp.gff");
    system("perl $scriptDir/compare-predictions-pdiff.pl $fileA cacpp-temp.gff $outDir $scriptDir");
    system("rm cacpp-temp.gff");
  }
}
