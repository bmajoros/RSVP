#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "cat file | $name fieldname > outfile\n" unless @ARGV==1;
my ($fieldName)=@ARGV;

while(<STDIN>) {
  chomp;
  my @fields=split/\s+/,$_;
  my $extraFields=pop @fields;
  $extraFields=~/$fieldName\=([^;]+)/;
  my $value=$1;
  print "$value\n";
}




