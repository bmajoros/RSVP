#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;
use Getopt::Std;

our $opt_s;
getopts('s');

my $name=ProgramName::get();
die "$name <in.gff> <out.bed>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($infile);

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
foreach my $transcript (@$transcripts) {
    # chrY    2770206 2770283 ENSG00000129824.6_2     0       +
    my $chrom=$transcript->getSubstrate();
    unless($chrom=~/chr/) {$chrom="chr$chrom"}
    my $strand=$transcript->getStrand();
    my $transcriptID=$transcript->getTranscriptId();
    my $numExons=$transcript->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i) {
        my $exon=$transcript->getIthExon($i);
        my $begin=$exon->getBegin();
        my $end=$exon->getEnd();
        if($opt_s) {
            my $exonType=$exon->getType();
            my $hasStart=
              ($exonType=~/initial-exon/ || $exonType=~/single-exon/);
            my $hasStop=
              ($exonType=~/final-exon/ || $exonType=~/single-exon/);
            if($hasStart) {
                if($strand eq "+") {$begin+=3}
                else {$end-=3}
            }
            if($hasStop) {
                if($strand eq "+") {$end-=3}
                else {$begin+=3}
            }
        }
        next unless $begin<$end;
        my $exonID="${transcriptID}_$i";
        print OUT "$chrom\t$begin\t$end\t$exonID\t0\t$strand\n";
    }
}
close(OUT);
