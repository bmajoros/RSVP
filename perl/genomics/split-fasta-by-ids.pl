#!/usr/bin/perl
use strict;

# Splits a multi-fasta file into multiple files, one file per sequence,
# using the sequence identifiers (first token after the '>') as filenames

my $usage="$0 <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $isOpen=0;
open(IN,$infile) || die "can't open $infile\n";
while(<IN>)
{
    if(/^\s*>\s*(\S+)/)
    {
	my $id=$1;
	if($isOpen) {close(OUT)}
	my $outfile="$id.fasta";
	open(OUT,">$outfile") || die "can't create file $outfile\n";
	$isOpen=1;
    }
    if($isOpen) {print OUT}
}
close(IN);
if($isOpen) {close(OUT)}


