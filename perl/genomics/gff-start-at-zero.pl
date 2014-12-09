#!/usr/bin/perl
use strict;
use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($filename);
my $n=@$transcripts;
open(OUT,">$filename") || die "can't write to $filename\n";
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    $transcript->shiftCoords(-$transcript->getBegin());
    my $gff=$transcript->toGff();
    print OUT $gff;
  }
close(OUT);


