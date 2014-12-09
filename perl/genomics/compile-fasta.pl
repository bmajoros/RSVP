#!/usr/bin/perl
use strict;

die "compile-fasta.pl <infile> <outfile>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't create $outfile\n";
open(IN,$infile) || die "can't open $infile\n";
my $defline=<IN>;
while(<IN>)
  {
    $_=~s/\s+//g;
    print OUT $_;
  }
close(IN);
close(OUT);


