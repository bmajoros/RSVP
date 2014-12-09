#!/usr/bin/perl
use strict;

# Swaps phases 1 and 2 so annotations appear correctly in IGV browser
# (due to apparent bug in IGV / nontraditional phase definition in TAIR)

while(<STDIN>) {
  chomp;
  my @fields=split/\t/,$_;
  next unless @fields>=8;
  if($fields[7]>0) {
    $fields[7]=3-$fields[7];
  }
  my $line=join("\t",@fields);
  print "$line\n";
}

