#!/usr/bin/perl
use strict;

use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($filename);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $n=$transcript->numExons();
    my $id=$transcript->getID();
    for(my $j=0 ; $j<$n ; ++$j)
      {
	my $exon=$transcript->getIthExon($j);
	my $length=$exon->getLength();
	print "$length : exon #$j in gene $id\n";
      }
  }





