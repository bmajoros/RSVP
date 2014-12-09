#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use SummaryStats;

$0=~/([^\/]+)$/;
my $usage="$1 <*.gff> <num-genes> <selected.gff> <not-selected.gff>";
die "$usage\n" unless @ARGV==4;
my ($filename,$N,$out1,$out2)=@ARGV;

# read the transcripts from the GFF file
my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($filename);
my $numTranscripts=@$transcripts;
if($N>$numTranscripts || $N<0) 
  {die "N must be between 0 and $numTranscripts\n"}

# randomly shuffle the transcripts
for(my $i=0 ; $i<$numTranscripts-1 ; ++$i)
  {
    my $range=$numTranscripts-$i;
    my $j=int(rand($range))+$i;
    my $tmp=$transcripts->[$i];
    $transcripts->[$i]=$transcripts->[$j];
    $transcripts->[$j]=$tmp;
  }

# write the first N into the first file
open(OUT1,">$out1") || die;
my $i;
for($i=0 ; $i<$N ; ++$i)
  {
    my $gff=$transcripts->[$i]->toGff();
    print OUT1 $gff;
  }
close(OUT1);

# write the rest into the second file
open(OUT2,">$out2") || die;
for(; $i<$numTranscripts ; ++$i)
  {
    my $gff=$transcripts->[$i]->toGff();
    print OUT2 $gff;
  }
close(OUT2);

