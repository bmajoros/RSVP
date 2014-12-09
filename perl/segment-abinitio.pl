#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use FastaWriter;
use TempFilename;
use GffTranscriptReader;
$|=1;

my $name=ProgramName::get();
die "$name <genezilla.iso> <chromosomes.fasta> <window-size> <padding>  <package-dir>\n"
  unless @ARGV==5;
my ($isoFile,$genomeFile,$windowSize,$padding,$packageDir)=@ARGV;

# Some initialization
my $tempFile="segment.tmp1";#TempFilename::generate();
my $tempFile2="segment.tmp2";#TempFilename::generate();
my $tempFile3="segment.tmp3";#TempFilename::generate();
my $genezilla="$packageDir/c++/genezilla/genezilla-abinitio";
die "Please run make-without-graphs\n" unless -e $genezilla;
my $writer=new FastaWriter;

# Read the chromosomes, process one by one
my $reader=new FastaReader($genomeFile);
my $gffReader=new GffTranscriptReader();
while(1) {
  my ($defline,$chrSeq)=$reader->nextSequence();
  last unless defined $defline;
  $chrSeq="\U$chrSeq";
  $chrSeq=~s/[^ACGT]//g;
  my $chrLen=length($chrSeq);
  $defline=~/^\s*>\s*(\S+)/ || die "Can't parse defline: $defline\n";
  my $chrID=$1;

  # Process each overlapping chunk, collect gene predictions
  my $increment=int($windowSize/2);
  my @genes;
  for(my $begin=0 ; $begin<$chrLen-10 ; $begin+=$increment) {
    my $end=$begin+$windowSize;
    if($end>$chrLen) { $end=$chrLen }
    my $chunkLen=$end-$begin;
    my $chunkSeq=substr($chrSeq,$begin,$chunkLen);
    my $chunkID="${chrID}:$begin";

    # Write chunk into file
    $writer->writeFastaFromRef(">$chunkID /chr=$chrID /offset=$begin",
			      \$chunkSeq,$tempFile);

    # Run genezilla on chunk
    my $cmd="$genezilla $isoFile $tempFile";
    system("bash -c '$cmd 1> $tempFile2 2> $tempFile3'");
    #system("bash -c '$cmd 1> $tempFile2 '");
    my $stderr=`cat $tempFile3`;
    die $stderr unless $stderr=~/TOTAL MEMORY USAGE/;
    my $transcripts=$gffReader->loadGFF($tempFile2);
    foreach my $t (@$transcripts) {$t->shiftCoords($begin)}
    push @genes,@$transcripts;
  }

  # Resolve overlapping gene predictions
  @genes=sort {$a->getBegin() <=> $b->getBegin()} @genes;
  my $numTranscripts=@genes;
  for(my $i=0 ; $i<$numTranscripts ; ) {
    my $j;
    for($j=$i ; $j<$numTranscripts-1 ; ++$j)
      { last unless $genes[$j]->overlaps($genes[$j+1]) }

    ###
    $j=$i;
    ###

    my $begin=$genes[$i]->getBegin()-$padding;
    my $end=$genes[$j]->getEnd()+$padding;
    if($begin<0) { $begin=0 }
    if($end>$chrLen) { $end=$chrLen }
    my $numToMerge=$j-$i+1;
    print "$chrID\tGZ\tchunk\t$begin\t$end\t$numToMerge\t+\t.\n";
    $i=$j+1;
  }
}
$reader->close();

system("rm -f $tempFile $tempFile2 $tempFile3");

