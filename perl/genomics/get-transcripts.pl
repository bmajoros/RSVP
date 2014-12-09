#!/usr/bin/perl
#################################################################
# get-transcripts.pl
# bmajoros@tigr.org, July 2002
# This software is govered by the Gnu Public License (GPL).
#################################################################
use strict;
use SummaryStats;
use Translation;
use FastaWriter;
use FastaReader;
use FileHandle;
use GffTranscriptReader;
use Transcript;
use Getopt::Std;
use ProgramName;
$|=1;

our($opt_e);
getopts("e");
my $perExon=$opt_e;

my $name=ProgramName::get();
my $usage="
$name [-e] <*.gff> <in.fasta> <out.fasta>

-e = separate exons into individual FASTA sequences

";
die "$usage\n" unless @ARGV==3;
my ($gffFilename,$seqFilename,$outFilename)=@ARGV;

# Read the sequences from the FASTA file
my %sequences;
my $fastaReader=new FastaReader($seqFilename);
while(1)
  {
    my ($defline,$seqRef)=$fastaReader->nextSequenceRef();
    last unless defined $defline;
    $defline=~/>\s*(\S+)/;
    my $id=$1;
    $sequences{$id}=$seqRef;
  }
$seqFilename=~/([^\/]+)$/;
my $seqBase=$1;

# Read the feature coordinates from the GFF file
my $gffReader=new GffTranscriptReader;
$gffReader->doNotSortTranscripts();
my $transcripts=$gffReader->loadGFF($gffFilename);

my $fastaWriter=new FastaWriter;
my $filehandle=new FileHandle(">$outFilename") ||
    die "Can't create $outFilename";

my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $trans=$transcripts->[$i];
    my $substrate=$trans->{substrate};
    my $genomicSequence=$sequences{$substrate} || die $substrate;
    my $transcript=$trans->loadTranscriptSeq($genomicSequence);
    my $transcriptId=$trans->{transcriptId};
    my $numExons=@{$trans->{exons}};
    my $strand=$trans->{strand};
    my $begin=$trans->{begin};
    my $end=$trans->{end};
    my $startCodon=$trans->{startCodon};
    my $startCodonLexeme=substr($transcript,$startCodon,3);###
    my $isPartial=($trans->isPartial() ? "true" : "false");
    if($perExon)
      {
	for(my $i=0 ; $i<$numExons ; ++$i)
	  {
	    my $exon=$trans->getIthExon($i);
	    my $seq=$exon->getSequence();
	    $fastaWriter->addToFasta(">$transcriptId /exonNum=$i /numExons=$numExons /strand=$strand /begin=$begin /end=$end /substrate=$substrate /seqFile=$seqBase /startCodon=$startCodon /coordsystem=\"0-based,space-based\" /partial=$isPartial",$seq,$filehandle);
	  }
      }
    else
      {
	$fastaWriter->addToFasta(">$transcriptId /numExons=$numExons /strand=$strand /begin=$begin /end=$end /seqFile=$seqBase /startCodon=$startCodon /coordsystem=\"0-based,space-based\" /partial=$isPartial",$transcript,$filehandle);
      }
    undef $trans->{sequence};
  }
close($filehandle);




