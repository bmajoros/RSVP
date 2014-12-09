#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <filename>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my %exonTypes;
$exonTypes{"Sngl"}="single-exon";
$exonTypes{"Init"}="initial-exon";
$exonTypes{"Intr"}="internal-exon";
$exonTypes{"Term"}="final-exon";

my $substrate;
open(IN,$infile) || die "can't open $infile\n";

print "##gff-version 2\n";
while(<IN>)
  {
    if(/^(GENSCAN\S*)\s+(\S+)/)
      {
	print "##source-version $1 $2\n";
      }
    if(/Date run:\s*(\S+)\s*Time:\s*(\S+)/)
      {
	print "##date $1 $2\n\n";
      }
    if(/^Sequence\s+(\S+)/)
      {
	$substrate=$1;
      }
    if(/^\s*(\d+)\.\d+\s+([A-Za-z]+)\s+(\S)\s+(\d+)\s+(\d+)\s+\d+\s+(\d+)\s+.*\s+(\S+)\s*$/)
      {
	my ($transId,$exonType,$strand,$begin,$end,$frame,$score)=
	  ($1,$2,$3,$4,$5,$6,$7);
	if($begin>$end) {($begin,$end)=($end,$begin)}
	if($exonTypes{$exonType})
	  {
	    $exonType=$exonTypes{$exonType};
	    print "$substrate\tgenscan\t$exonType\t$begin\t$end\t$score\t$strand\t$frame\ttransgrp=$transId;\n";
	  }
      }
  }
close(IN);









