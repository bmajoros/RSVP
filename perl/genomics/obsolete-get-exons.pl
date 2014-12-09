#!/usr/bin/perl
use strict;
use lib('/home/bmajoros/perlib');
use SummaryStats;
use Translation;
use FastaWriter;
use FastaReader;
use FileHandle;
$|=1;

my $usage="$0 <*.gff> <*.seq>";
die "$usage\n" unless @ARGV==2;
my ($gffFilename,$seqFilename)=@ARGV;

# Read the sequence from the FASTA file
my $reader=new FastaReader($seqFilename);
my ($defline,$sequence)=$reader->nextSequence();

$seqFilename=~/([^\/]+)$/;
my $seqBase=$1;

# Process the GFF file
my %transcripts;
open(GFF,$gffFilename) || die $gffFilename;
while(<GFF>)
  {
    if(/exon\s+(\d+)\s+(\d+).*([\-\+]).*transgrp=(\S+)/)
      {
	my ($exonBegin,$exonEnd,$strand,$transcriptId)=($1-1,$2,$3,$4);
	if($exonBegin>$exonEnd)
	  {($exonBegin,$exonEnd)=($exonEnd,$exonBegin)}
	my $transcript=$transcripts{$transcriptId};
	if(!defined $transcript) 
	  {
	    $transcripts{$transcriptId}=$transcript=
	      {
	       transcriptId=>$transcriptId,
	       strand=>$strand
	      };
	  }
	my $exon=[$exonBegin,$exonEnd];
	if($strand eq "+")
	  {
	    push @{$transcript->{exons}},$exon;
	    if(!defined $transcript->{fivePrime})
	      {$transcript->{fivePrime}=$exonBegin}
	    $transcript->{threePrime}=$exonEnd;
	  }
	else
	  {
	    die unless $strand eq "-";
	    unshift @{$transcript->{exons}},$exon;
	    if(!defined $transcript->{threePrime})
	      {$transcript->{threePrime}=$exonEnd}
	    $transcript->{fivePrime}=$exonBegin;
	  }
      }
  }
close(GFF);

my (%donors,%acceptors,$numIntrons);
my @transcripts=values %transcripts;
undef %transcripts;
my $numTranscripts=@transcripts;
my $uniqueId=1000;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    ++$uniqueId;
    my $trans=$transcripts[$i];
    my $exons=$trans->{exons};
    my $transcriptId=$trans->{transcriptId};
    my $numExons=@$exons;
    my $strand=$trans->{strand};
    
    my $transcript;
    my @displayList;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	# Get next exon:
	my $exon=$exons->[$i];
	my $begin=$exon->[0];
	my $end=$exon->[1];
	my $exonLen=$end-$begin;
	my $exonSeq=substr($sequence,$begin,$exonLen);
	if($strand eq "-")
	  {
	    $exonSeq=Translation::reverseComplement(\$exonSeq);
	    unshift @displayList,$exonSeq;
	  }
	else {push @displayList,$exonSeq}

	# Get next intron:
	if($i<$numExons-1)
	  {
	    my $intronBegin=$end;
	    my $nextExon=$exons->[$i+1];
	    my $intronEnd=$nextExon->[0];
	    my $intronLen=$intronEnd-$intronBegin;
	    my $intron=substr($sequence,$intronBegin,$intronLen);
	    if($strand eq "-")
	      {
		$intron=Translation::reverseComplement(\$intron);
		unshift @displayList,$intron;
	      }
	    else {push @displayList,$intron}
	  }
      }

    #print "="x60 . "\n";
    #print "TRANSCRIPT: $transcriptId ($strand)\n";
    for(my $i=0 ; $i<@displayList ; $i+=2)
      {
	my $exonSeq=$displayList[$i];
	print "$exonSeq\n";
	next unless $i<$numExons-1;
	my $intronSeq=$displayList[$i+1];
	#print "INTRON $i: $intronSeq\n";# unless $intronSeq=~/^GT.*AG$/;
	$intronSeq=~/^(..).*(..)$/;
	my ($donor,$acceptor)=($1,$2);
	++$donors{$donor};
	++$acceptors{$acceptor};
	++$numIntrons;
      }
  }

exit;

my @donors=keys %donors;
my @acceptors=keys %acceptors;
foreach my $donor (@donors)
  {
    my $n=$donors{$donor};
    my $p=int(100*$n/$numIntrons+0.5);
    print "DONOR: $donor ($p\%)\n";
  }
foreach my $acceptor (@acceptors)
  {
    my $n=$acceptors{$acceptor};
    my $p=int(100*$n/$numIntrons+0.5);
    print "ACCEPTOR: $acceptor ($p\%)\n";
  }




