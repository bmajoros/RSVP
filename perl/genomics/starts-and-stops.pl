#!/usr/bin/perl
use strict;
use FastaReader;

die "starts-and-stops.pl <transcripts.fasta> <starts.out> <stops.out>\n" 
  unless @ARGV==3;
my ($infile,$starts,$stops)=@ARGV;

open(STARTS,">$starts") || die "Can't write to file: $starts\n";
open(STOPS,">$stops") || die "Can't write to file: $stops\n";
my $reader=new FastaReader($infile);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    $def=~/^>(\S+)/;
    my $id=$1;
    my $start=substr($seq,0,3);
    my $stop=substr($seq,length($seq)-3,3);
    print STARTS "$start $id\n";
    print STOPS "$stop $id\n";
  }
close(STOPS);
close(STARTS);


