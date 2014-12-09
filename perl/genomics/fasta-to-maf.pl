#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <in.fasta> <out.maf>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "Can't write into file: $outfile\n";
print OUT "##maf version=12\n
a score=1\n";

my $reader=new FastaReader($infile);
while(1)
{
    my ($defline,$seq)=$reader->nextSequenceRef();
    last unless defined $defline;
    $defline=~/^>(\S+)/ 
      || die "can't parse defline: $defline\n";
    my $species=$1;
    $species=~s/_aligned$//g;
    my $length=length($$seq);
    print OUT "s $species 0 $length + $length $$seq\n";
}
print OUT "\n";
close(OUT);
