#!/usr/bin/perl
use strict;
use GffReader;
use TempFilename;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.gff> <new-substrate-id>";
die "$usage\n" unless @ARGV==2;
my ($infile,$substrate)=@ARGV;

my $tempfile=TempFilename::generate();
open(OUT,">$tempfile") || die "Can't create temp file: $tempfile\n";
my $reader=new GffReader;
my $features=$reader->loadGFF($infile);
my $numFeatures=@$features;
for(my $i=0 ; $i<$numFeatures ; ++$i)
  {
    my $feature=$features->[$i];
    $feature->setSubstrate($substrate);
    my $gff=$feature->toGff();
    print OUT $gff;
  }
close(OUT);
rename($tempfile,$infile) || die "Can't rename $tempfile as $infile\n";






