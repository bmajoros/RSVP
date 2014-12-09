#!/usr/bin/perl
#Given a folder of GFF files, retain only lines with fullSupport "1" or fullSupport=1
#and output the resulting files to another folder.

use warnings;
use strict;

die "keep-supported-predictions.pl <gff-dir> <out-dir>" unless @ARGV==2;

my ($gffDir, $outDir) = @ARGV;

my $gffs = `ls $gffDir/*.gff`;
my @gffList = split(/\s+/, $gffs);
for my $gff (@gffList)
{
  my ($gffName) = ($gff =~ /(\w+).gff/);
  system("awk '/fullSupport \"1\"/ || /fullSupport=1/ {print}' $gffDir/$gffName.gff > $outDir/$gffName.gff");
}
