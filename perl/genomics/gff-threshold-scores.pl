#!/usr/bin/perl
use strict;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.gff> <score>\n" unless @ARGV==2;
my ($infile,$minScore)=@ARGV;

my $reader=new GffReader;
my $features=$reader->loadGFF($infile);
my $n=@$features;
for(my $i=0 ; $i<$n ; ++$i) {
  my $feature=$features->[$i];
  if($feature->getScore()>=$minScore) { print $feature->toGff() }
}


