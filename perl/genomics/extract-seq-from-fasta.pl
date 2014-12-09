#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

$0=~/([^\/]+)$/;
my $usage="$1 <seq-number> <filename>";
die "$usage\n" unless @ARGV==2;
my ($index,$infile)=@ARGV;

my $reader=new FastaReader($infile);
for(my $i=0 ; $i<$index ; ++$i)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
  }

my $writer=new FastaWriter;

my ($defline,$sequence)=$reader->nextSequence();
if(!defined($defline)) {die "No such sequence in fasta file.\n"}
if($defline=~/\n/) {chop $defline}
if($sequence=~/\n/) {chop $sequence}
$writer->addToFasta($defline,$sequence,\*STDOUT);


