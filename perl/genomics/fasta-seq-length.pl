#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $len=0;
open(IN,$filename) || die;
while(<IN>)
  {
    unless(/>/)
      {
	#$len+=($_=~s/([ATCGNatcgn])/$1/g);
	$len+=($_=~s/([A-Za-z])/$1/g);
      }
  }
close(IN);
print "$len bp\n";


