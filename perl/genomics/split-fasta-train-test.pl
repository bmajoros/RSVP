#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <infile.fasta> <num-train> <out-train.fasta> <out-test.fasta>";
die "$usage\n" unless @ARGV==4;
my ($infile,$numTrain,$trainFile,$testFile)=@ARGV;

my $reader=new FastaReader($infile);
my $writer=new FastaWriter;

open(TRAIN,">$trainFile") || die "can't create $trainFile\n";
for(my $i=0 ; $i<$numTrain ; ++$i)
{
    my ($defline,$seq)=$reader->nextSequence();
    $writer->addToFasta($defline,$seq,\*TRAIN);
}
close(TRAIN);

open(TEST,">$testFile") || die "can't create $testFile\n";
while(1)
{
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined $defline;
    $writer->addToFasta($defline,$seq,\*TEST);
}
close(TEST);
