#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "cat file | $name fieldname > outfile\n" unless @ARGV==1;
my ($fieldName)=@ARGV;
my @fields=("substrate","source","type","begin","end","score","strand","frame","extra");
my $index;
for(my $i=0 ; $i<@fields ; ++$i) { if($fieldName eq $fields[$i]) {$index=$i} }
if(!defined($index)) {
  die "field must be one of: substrate|source|type|begin|end|score|strand|frame|extra\n";
}

while(<STDIN>) {
  chomp;
  if(/^\s*#/) {next}
  my @fields=split/\s+/,$_;
  my $value=$fields[$index];
  print "$value\n";
}



