#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile> \n" unless @ARGV==1;
my ($infile)=@ARGV;

my %seen;
my $reader=new FastaReader($infile);
my $writer=new FastaWriter();
while(1) {
  my ($defline,$sequence)=$reader->nextSequence();
  last unless $defline;
  next if $seen{$sequence};
  $seen{$sequence}=1;
  $writer->addToFasta($defline,$sequence,\*STDOUT);
}
$reader->close();

