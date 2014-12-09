#!/usr/bin/perl
use strict;
use ProgramName;
use GffTranscriptReader;

my $name=ProgramName::get();
my $usage="$name <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($infile);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $numExons=$transcript->numExons();
    my $lastExon=$transcript->getIthExon($numExons-1);
    $lastExon->trimFinalPortion(3);
    my $gff=$transcript->toGff();
    print $gff;
  }



