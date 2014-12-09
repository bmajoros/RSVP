#!/usr/bin/perl
use strict;

my $usage="$0 <infile.gff>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
while(<IN>)
  {
    next unless $_=~/\S/;
    my @fields=split/\s+/,$_;
    ++$fields[3];
    my $line=join("\t",@fields);
    print "$line\n";
  }
close(IN);

