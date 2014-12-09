#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "cat file | $name > outfile\n" unless @ARGV==0;

while(<STDIN>) {
  chomp;
  if(/^\s*#/) {next}
  my @fields=split/\s+/,$_;
  my $length=$fields[4]-$fields[3];
  print "$length\n";
}



