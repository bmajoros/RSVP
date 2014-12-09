#!/usr/bin/perl
##################################################################
#
# Reads a FASTA file from STDIN and outputs it to STDOUT with
# proper formatting (60 bases wide)
#
##################################################################
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use FastaReader;
use FastaWriter;


my $usage="cat infile.fasta | format-fasta.pl N/P > outfile.fasta";
die "$usage\n" unless @ARGV==1;
my ($type)=@ARGV;

my $writer=new FastaWriter;
my $reader=FastaReader::readerFromFileHandle(\*STDIN);
while(1)
  {
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined $defline;
    if($type eq "N") {$seq=~s/[^ATCGN]/N/g}
    else {$seq=~s/[^ARNDCQEGHILKMFPSTWYV]/N/g}
    $writer->addToFasta($defline,$seq,\*STDOUT);
  }




