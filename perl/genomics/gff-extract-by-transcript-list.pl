#!/usr/bin/perl
use strict;
use GffTranscriptReader;

# Process command line and load GFF file
$0=~/([^\/]+)\s*$/;
my $usage="$1 <in.gff> <genelist.txt> <out.gff>";
die "$usage\n" unless @ARGV==3;
my ($inGff,$geneFile,$outGff)=@ARGV;
my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($inGff);

# Load gene ID's into hash table
my %keep;
open(IN,$geneFile) || die "can't open $geneFile\n";
while(<IN>)
  {
    if(/(\S+)/) {$keep{$1}=1}
  }
close(IN);

# Output the transcripts having a gene ID in the list
open(OUT,">$outGff") || die "can't write to $outGff\n";
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $id=$transcript->getTranscriptId();
    next unless $keep{$id};
    my $gff=$transcript->toGff();
    print OUT $gff;
  }
close(OUT);
