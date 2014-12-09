#!/usr/bin/perl
use strict;
use FastaReader;

$0=~/([^\/]+)\s*$/;
die "$1 <*.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    print "$seq\n";
  }
