#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

$0=~/([^\/]+)$/;
my $usage="$1 <infile> <outfile> <N>";
die "$usage\n" unless @ARGV==3;
my ($infile,$outfile,$n)=@ARGV;

open(OUT,">$outfile") || die "can't create $outfile\n";
my $writer=new FastaWriter;
my $reader=new FastaReader($infile);
my $count=0;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    $writer->addToFasta($def,$seq,\*OUT);
    ++$count;
  }

close(OUT);
print "$count sequences extracted\n";



