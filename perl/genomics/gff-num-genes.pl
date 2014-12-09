#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <*.gff>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGenes($infile);
my $n=@$transcripts;
print "$n\n";

