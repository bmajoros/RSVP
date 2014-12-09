#!/usr/bin/perl
use strict;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
my $usage="cat seq.txt | $name <identifier>";
die "$usage\n" unless @ARGV==1;
my ($id)=@ARGV;

my $seq;
while(<STDIN>)
  {
    $_=~s/\s+//g;
    $seq.=$_;
  }

my $writer=new FastaWriter;
$writer->addToFasta(">$id",$seq,\*STDOUT);

