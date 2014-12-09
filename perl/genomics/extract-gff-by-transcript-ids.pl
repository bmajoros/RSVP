#!/usr/bin/perl
use strict;
use GffTranscriptReader;

my $usage="$0 <*.gff> <id-file>";
die "$usage\n" unless @ARGV==2;
my ($gffFilename,$idFilename)=@ARGV;

my %keep;
open(IN,$idFilename) || die "can't open $idFilename";
while(<IN>)
  {
    chop;
    if(/(\S+)/)
      {
	#print "==============>[$1]\n";
	$keep{$1}=1;
      }
  }
close(IN);

my $transcriptReader=new GffTranscriptReader;
my $transcripts=$transcriptReader->loadGFF($gffFilename);

foreach my $transcript (@$transcripts)
  {
    my $id=$transcript->getID();
    $id=~/([^;]+)/;
    $id=$1;
    #print "[$id]\n";
    next unless $keep{$id};
    my $gff=$transcript->toGff();
    print "$gff\n";
  }



