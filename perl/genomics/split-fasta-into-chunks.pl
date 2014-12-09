#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

$0=~/([^\/]+$)/;
die "$1 <infile> <outfile> <chunk-size>\n" unless @ARGV==3;
my ($infile,$outfile,$chunkSize)=@ARGV;

open(OUT,">$outfile") || die "Can't create file: $outfile\n";
my $writer=new FastaWriter();
my $reader=new FastaReader($infile);
while(1)
  {
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined $defline;
    $defline=~/^\s*>\s*([a-zA-Z01-9_-]+)(.*)/ || 
      die "can't parse defline: $defline\n";
    my ($id,$rest)=($1,$2);
    my $seqLen=length $seq;
    my $numChunks=$seqLen/$chunkSize;
    my $pos=0;
    for(my $i=0 ; $i<$numChunks ; ++$i)
      {
	my $thisChunkSize=($chunkSize<$seqLen ? $chunkSize : $seqLen);
	my $chunk=substr($seq,$pos,$thisChunkSize);
	$seqLen-=$thisChunkSize;
	$pos+=$thisChunkSize;
	my $chunkId=$i+1;
	my $newDef=">$id-$chunkId$rest";
	$writer->addToFasta($newDef,$chunk,\*OUT);
      }
    print STDOUT "$numChunks chunks written\n";
  }
close(OUT);


