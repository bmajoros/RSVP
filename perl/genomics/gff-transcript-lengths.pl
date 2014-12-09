#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($filename);
my $n=@$transcripts;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $length=$transcript->getExtent();
    print "$length\n";
  }





