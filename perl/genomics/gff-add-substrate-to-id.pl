#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;
use TempFilename;

my $name=ProgramName::get();
my $usage="$name <in-and-out.gff>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $tmp=TempFilename::generate();
open(OUT,">$tmp") || die "Can't create temp file: $tmp\n";
my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($infile);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $id=$transcript->getID();
    my $substrate=$transcript->getSubstrate();
    $id="$substrate-$id";
    $transcript->setTranscriptId($id);
    $transcript->setGeneId($id);
    my $gff=$transcript->toGff();
    print OUT "$gff\n";
  }
close(OUT);
system("mv $tmp $infile");




