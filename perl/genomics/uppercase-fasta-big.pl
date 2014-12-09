#!/usr/bin/perl
use strict;
use TempFilename;

my $usage="$0 <filename>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $tmpFile=TempFilename::generate();
open(OUT,">$tmpFile");
open(IN,$infile) || die "can't open $infile\n";
while(<IN>)
  {
    if(/^\s*>/) {print OUT}
    else
      {
	$_="\U$_";
	print OUT;
      }
  }
close(IN);
close(OUT);

system("mv $tmpFile $infile");

