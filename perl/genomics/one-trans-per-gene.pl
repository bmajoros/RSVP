#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> <out.gff>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $reader=new GffTranscriptReader();
my $geneList=$reader->loadGenes($infile);

open(OUT,">$outfile") || die "can't write to file: $outfile\n";
foreach my $gene (@$geneList) {
    my $transcript=$gene->getIthTranscript(0);
    my $gff=$transcript->toGff();
    print OUT "$gff\n";
}
close(OUT);









