#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;
$|=1;

my $name=ProgramName::get();
die "$name <segments.gff> <pileup-dir> <package-dir>\n" unless @ARGV==3;
my ($infile,$pileupDir,$packageDir)=@ARGV;

my $reader=new GffReader;
my $bySubstrateHash=$reader->hashBySubstrate($infile);
my @substrates=keys %$bySubstrateHash;
foreach my $substrate (@substrates) {
  print STDERR "processing $substrate...\n";
  my $cmd="$packageDir/c++/dump-pileup $pileupDir/$substrate.pileup";
  my @pileup=`$cmd`;
  my $chrLen=@pileup;
  my $segments=$bySubstrateHash->{$substrate};
  my $numSegments=@$segments;
  for(my $i=0 ; $i<$numSegments ; ++$i) {
    my $segment=$segments->[$i];
    my $begin=$segment->getBegin();
    my $end=$segment->getEnd();
    while($begin>=0 && $pileup[$begin]>0) { --$begin }
    while($end<$chrLen && $pileup[$end]>0) { ++$end }
    if($begin<$segment->getBegin()) {
      my $extend=$segment->getBegin()-$begin;
      print "extended begin by $extend\n";
    }
    if($segment->getEnd()<$end) {
      my $extend=$end-$segment->getEnd();
      print "extended end by $extend\n";
    }
  }
}



