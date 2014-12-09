#!/usr/bin/perl
use strict;
use GffTranscriptReader;
$|=1;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;
use SummaryStats;

my $reader=new GffTranscriptReader;
my $transcripts=$reader->loadGFF($filename);
my $n=@$transcripts;
my (@transcriptLengths,@transcriptExtents,@exonCounts,$numPlusStrand,
    @transcriptScores);
foreach my $transcript (@$transcripts)
  {
    push @transcriptLengths,$transcript->getLength();
    push @transcriptExtents,$transcript->getExtent();
    push @exonCounts,$transcript->numExons();
    if($transcript->getStrand() eq "+"){++$numPlusStrand}
    push @transcriptScores,$transcript->getScore();
  }
my $percentPlusStrand=int(1000*$numPlusStrand/$n+5/9)/10;

my ($meanLen,$stddevLen,$minLen,$maxLen)=
  SummaryStats::roundedSummaryStats(\@transcriptLengths);
my ($meanExt,$stddevExt,$minExt,$maxExt)=
  SummaryStats::roundedSummaryStats(\@transcriptExtents);
my ($meanExn,$stddevExn,$minExn,$maxExn)=
  SummaryStats::roundedSummaryStats(\@exonCounts);
my ($meanScr,$stddevScr,$minScr,$maxScr)=
  SummaryStats::roundedSummaryStats(\@transcriptScores);

print "$n transcripts, $percentPlusStrand% on + strand\n";
print "transcript length: $meanLen+/-$stddevLen ($minLen-$maxLen)\n";
print "transcript extent (including introns): $meanExt+/-$stddevExt ($minExt-$maxExt)\n";
print "exons per transcript: $meanExn+/-$stddevExn ($minExn-$maxExn)\n";
print "transcript score: $meanScr+/-$stddevScr ($minScr-$maxScr)\n";
