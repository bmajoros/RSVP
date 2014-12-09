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
print "$numTranscripts transcripts\n";

my %substrates;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    my $substrate=$transcript->getSubstrate();
    push @{$substrates{$substrate}},$transcript;
  }

my @substrates=keys %substrates;
foreach my $substrate (@substrates)
  {
    my $transcripts=$substrates{$substrate};
    my $outfile="$substrate.gff";
    open(OUT,">$outfile") || die "Can't creat $outfile\n";
    foreach my $transcript (@$transcripts)
      {
	print OUT $transcript->toGff();
      }
    close(OUT);
  }
