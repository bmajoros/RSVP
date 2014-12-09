#!/usr/bin/perl
use strict;

my $usage="$0 <infile> <contig-ID>";
die "$usage\n" unless @ARGV==2;
my ($infile,$substrate)=@ARGV;

my %types;
$types{"Internal"}="internal-exon";
$types{"Initial"}="initial-exon";
$types{"Terminal"}="final-exon";
$types{"Single"}="single-exon";

open(IN,$infile) || die "can't open $infile\n";
while(<IN>)
  {
    if(/^\s*(\d+)\s+(\d+)\s*([+-])\s*([a-zA-Z]+)\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$/)
      {
	my ($geneID,$exonID,$strand,$type,$begin,$end)=
	  ($1,$2,$3,$4,$5,$6);
	my $exonType=$types{$type};
	die "unknown exon type: $type" unless defined $exonType;
	print "$substrate\tgenemark\t$exonType\t$begin\t$end\t.\t$strand\t.\ttransgrp=$geneID;\n";
      }
  }
close(IN);

