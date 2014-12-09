#!/usr/bin/perl
use strict;
use lib('/home/bmajoros/xml/xmlephant');
use xmlephant;

my $usage="$0 <xml-file>";
die "$usage\n" unless @ARGV==1;
my ($xmlFilename)=@ARGV;

my $query=xmlephant::parseQuery("/gene",$xmlFilename);
while(my $gene=$query->nextHit())
  {
    my $CT=$gene->lookupAttribute("CT");
    my $arm=$gene->lookupAttribute("arm");
    my $exons=$gene->findTags("exon");
    foreach my $exon (@$exons)
      {
	my $start=$exon->lookupAttribute("start");
	my $end=$exon->lookupAttribute("end");
	my $strand=($start<$end ? "+" : "-");
	print "$arm\tcelera\texon\t$start\t$end\t.\t$strand\t.\ttranscript=$CT\n";
      }
  }
