#!/usr/bin/perl
use strict;
use GenbankParser;
use FastaWriter;

$0=~/([^\/]+$)/;
my $program=$1;
my $usage="$program <infile>\n";
die $usage unless @ARGV==1;
my ($infile)=@ARGV;

my $transgrp=1;
my $parser=new GenbankParser($infile);
while(1)
  {
    my $entry=$parser->nextEntry();
    last unless $entry;

    my $defline=">$transgrp";
    my $organism=$entry->findUnique("ORGANISM");
      # || die "No ORGANISM field in Genbank entry!\n";
    if(!defined($organism)) {$organism="unknown_organism"}
    my @words=split/\s+/,$organism;
    $organism="$words[0]_$words[1]";
    my $seq=$entry->getSubstrate();

    my $geneId=$transgrp;
    my $locus=$entry->findUnique("LOCUS");
    if($locus) {$geneId=$locus}

    my $features=$entry->findUnique("FEATURES");
    if($features)
      {
	print "FEATURES:\n";
	my $n=$features->numPairs();
	for(my $i=0 ; $i<$n ; ++$i)
	  {
	    my $pair=$features->getIthPair($i);
	    my ($key,$value)=@$pair;
	    print "\t$key\n";
	  }
      }
  }




