#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcriptArray=$reader->loadGFF($infile);
my $numTranscripts=@$transcriptArray;

for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    my $id=$transcript->getID();
    $id=~s/;//g;
    my $outfile="$id.gff";
    open(OUT,">$outfile") || die "Can't create $outfile";
    print OUT $transcript->toGff();
    close(OUT);
  }
