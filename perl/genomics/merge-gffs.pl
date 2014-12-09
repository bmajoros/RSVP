#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;
use GffTranscriptReader;

my $name=ProgramName::get();
my $usage="$name <gff-dir> <map-file> <out.gff>";
die "$usage\n" unless @ARGV==3;
my ($dir,$mapFile,$outFile)=@ARGV;

# Load the mapping file
my %map;
open(IN,$mapFile) || die "Can't open $mapFile\n";
<IN>; # discard header line
while(<IN>)
  {
    if(/(\S+)\s+(\S+)\s+(\S+)/)
      {
	my ($chunkId,$contigId,$offset)=($1,$2,$3);
	$map{$chunkId}=[$contigId,$offset];
      }
  }
close(IN);

# Iterate through the GFF files in the directory
my $gffReader=new GffReader;
my $transcriptReader=new GffTranscriptReader;
my $ls=`ls $dir/*.gff`;
my @files=split/\s+/,$ls;
my $n=@files;
open(OUT,">$outFile") || die "can't create file: $outFile\n";
my $transcriptId=1;
my $featuresAreTranscripts=1;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $file=$files[$i];
    my $features=$transcriptReader->loadGFF($file);
    if(@$features==0)
      {
	$features=$gffReader->loadGFF($file);
	$featuresAreTranscripts=0;
      }
    my $m=@$features;
    for(my $j=0 ; $j<$m ; ++$j)
      {
	my $feature=$features->[$j];
	my $chunkId=$feature->getSubstrate();
	my $record=$map{$chunkId};
	if(!defined($record))
	  {
	    my $debug=$feature->toGff();
	    print "$debug\n";
	    print "chunkId=$chunkId is not defined in the map!";
	    die;
	  }
	my ($contigId,$offset)=@$record;
	$feature->shiftCoords($offset);
	$feature->setSubstrate($contigId);
	if($featuresAreTranscripts)
	  {
	    $feature->setTranscriptId($transcriptId);
	    $feature->setGeneId($transcriptId);
	    ++$transcriptId;
	  }
	my $gff=$feature->toGff();
	print OUT $gff;
      }
  }
close(OUT);

