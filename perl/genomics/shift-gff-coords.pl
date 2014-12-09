#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;

my $name=ProgramName::get();
die "$name <*.gff> <amount>\n" unless @ARGV==2;
my ($infile,$delta)=@ARGV;

my $reader=new GffReader();
my $featureArray=$reader->loadGFF($infile);
foreach my $feature (@$featureArray) {
  $feature->shiftCoords($delta);
  my $gff=$feature->toGff();
  print "$gff";
}



