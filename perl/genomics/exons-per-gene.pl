#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use SummaryStats;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($filename);

my @array;
foreach my $transcript (@$transcripts)
  {
    my $n=$transcript->numExons();
    push @array,$n;
    print "$n\n";
  }
my ($mean,$sd,$min,$max)=SummaryStats::roundedSummaryStats(\@array);
print STDERR "$mean +/- $sd ($min-$max)\n";



