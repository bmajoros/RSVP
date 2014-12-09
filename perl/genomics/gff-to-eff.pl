#!/usr/bin/perl
use strict;
use GffTranscriptReader;

die "$0 <infile>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($infile);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $strand=$transcript->getStrand();
    my $substrate=$transcript->getSubstrate();
    my $numExons=$transcript->numExons();
    for(my $j=0 ; $j<$numExons ; ++$j)
      {
	my $exon=$transcript->getIthExon($j);
	my $begin=$exon->getBegin();
	my $end=$exon->getEnd();
	if($strand eq "-") {($begin,$end)=($end,$begin)}
	print "$substrate $begin $end\n";
      }
    print "\n";
  }

