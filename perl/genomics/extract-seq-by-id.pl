#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl');
use FastaReader;
use FastaWriter;

my $usage="$0 <seq-id> <filename>";
die "$usage\n" unless @ARGV==2;
my ($id,$infile)=@ARGV;

my $writer=new FastaWriter;
my $reader=new FastaReader($infile);
my ($defline,$sequence);
while(1)
  {
    ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    $defline=~/^\s*>\s*([^\s;\|]+)/ || die;
    last if $1 eq $id;
  }

if(!defined($defline)) {die "$id: No such sequence in fasta file ($infile).\n"}
if($defline=~/\n/) {chop $defline}
if($sequence=~/\n/) {chop $sequence}
#print "$defline\n$sequence\n";
$writer->addToFasta($defline,$sequence,\*STDOUT);


