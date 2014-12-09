#!/usr/bin/perl
use strict;

#Run compare-predictions-asta.pl on all chunks in a folder

die "compare-all-chunk-predictions.pl <annotation-dir> <prediction-dir> <astalavista-path> <script-dir>" unless @ARGV==4;

print("pred.#\t\tann.#\t\tAStalavista\tEnds\n");

my ($dataDir, $outDir, $astaPath, $scriptDir) = @ARGV;

my $files = `ls $dataDir`;
my @files = split(/\s+/, $files);

for my $file (@files)
{
  if ($file =~ /(\d+).gff/)
  {
    my $fileA = "$dataDir/$1.gff";
    my $fileB = "$outDir/$1.gff";
    system("perl $scriptDir/keep-center-transcripts.pl $fileB > cacp-temp.gff");
    system("perl $scriptDir/compare-predictions-asta.pl $fileA cacp-temp.gff $astaPath");
    system("rm cacp-temp.gff");
  }
}
