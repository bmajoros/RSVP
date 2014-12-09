#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;
use Getopt::Std;

our $opt_s;
getopts('s');

my $name=ProgramName::get();
die "$name [-s] <in.gff> <out.bed> <max-length>
    where -s = omit splice sites\n" 
unless @ARGV==3;
my ($infile,$outfile,$maxLen)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($infile);

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
foreach my $transcript (@$transcripts) {
    my $chrom=$transcript->getSubstrate();
    unless($chrom=~/[cC]hr/) {$chrom="chr$chrom"}
    my $strand=$transcript->getStrand();
    my $transcriptID=$transcript->getTranscriptId();
    my $introns=$transcript->getIntrons();
    my $numIntrons=@$introns;
    for(my $i=0 ; $i<$numIntrons ; ++$i) {
        my $intron=$introns->[$i];
        my $begin=$intron->[0];
        my $end=$intron->[1];
        if($opt_s) {$begin+=2; $end-=2}
        if($end-$begin>$maxLen) {$end=$begin+$maxLen}
        my $intronID="${transcriptID}_${i}i";
        next unless $begin<$end;
        print OUT "$chrom\t$begin\t$end\t$intronID\t0\t$strand\n";
    }
}
close(OUT);
