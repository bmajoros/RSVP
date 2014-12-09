#!/usr/bin/perl
#Quick and dirty script to make a list of gene names from a chunk folder

use warnings;
use strict;

die "make-gene-name-list.pl <chunk-folder>" unless @ARGV==1;

my $chunkDir = shift;

my $list = `ls $chunkDir/*.coords.gff`;
my @files = split(/\s+/, $list);

for my $file (@files)
{
  open(FILE, $file);
  my $line = <FILE>;
  my @lineSplit = split(/\s+/, $line);
  my $geneString = $lineSplit[8];
  my $gene;
  if ($geneString =~ /gene_id=(\w+);/)
  {
    $gene = $1;
  }
  else
  {
    die "Error in finding gene from gene id";
  }
  print("$gene\n");
}
