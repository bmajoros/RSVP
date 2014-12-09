#!/usr/bin/perl
#################################################################
# gff-include-stop-codon.pl
# bmajoros@tigr.org, Feb 2004
# This software is govered by the Artistic License.
#################################################################
use strict;
use SummaryStats;
use Translation;
use FastaReader;
use FileHandle;
use GffTranscriptReader;
use Transcript;
$|=1;

use Getopt::Std;
our $opt_f;
our $opt_b;
getopts('fb');

$0=~/([^\/]+)\s*$/;
my $usage="$1 [options] <in.gff> <*.fasta> <out.gff> TAG,TGA,TAA
where -f=force final exon extension, even if no stop codon found
      -b=blind (fasta file is ignored; always adds 3 bases for the stop)";
die "$usage\n" unless @ARGV==4;
my ($gffFilename,$seqFilename,$outFilename,$stops)=@ARGV;

my @stops=split/[\s,\|]+/,$stops;
my %isStop;
foreach my $stop (@stops) {$isStop{$stop}=1}


# Read the sequences from the FASTA file
my %sequences;
if(!$opt_b) {
    my $fastaReader=new FastaReader($seqFilename);
    while(1) {
        my ($defline,$seqRef)=$fastaReader->nextSequenceRef();
        last unless defined $defline;
        $defline=~/>\s*(\S+)/;
        my $id=$1;
        $sequences{$id}=$seqRef;
    }
}
$seqFilename=~/([^\/]+)$/;
my $seqBase=$1;

# Read the feature coordinates from the GFF file
my $gffReader=new GffTranscriptReader;
$gffReader->setStopCodons(\%isStop);
$gffReader->doNotSortTranscripts();
my $transcripts=$gffReader->loadGFF($gffFilename);

my $gff;
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $trans=$transcripts->[$i];
    my $substrate=$trans->{substrate};
    my $genomicSequence=$opt_b ? undef : $sequences{$substrate} || die $substrate;
    my $substrateLen=$opt_b ? 0 : length $$genomicSequence;
    my $transSeq=$opt_b ? undef : $trans->loadTranscriptSeq($genomicSequence);
    my $transcriptId=$trans->{transcriptId};
    my $numExons=@{$trans->{exons}};
    my $strand=$trans->{strand};
    my $begin=$trans->{begin};
    my $end=$trans->{end};
    my $stopLexeme=$opt_b ? "XXX" : substr($transSeq,length($transSeq)-3,3);
    my $lastExon=$trans->getIthExon($numExons-1);
    my $lastExonType=$lastExon->getType();
    my $firstExon=$trans->getIthExon(0);
    my $firstExonType=$firstExon->getType();
    my $firstCodon=$opt_b ? "ATG" : substr($transSeq,0,3);
    if(($firstExonType eq "initial-exon" || $firstExonType eq "single-exon")
       && $firstCodon ne "ATG")
      {
	print "WARNING! $firstCodon is not ATG...deleting transcript\n";
	undef $transcripts->[$i];
	next;
      }
    if(($lastExonType eq "final-exon" || $lastExonType eq "single-exon") 
       && (!$isStop{$stopLexeme} || $opt_f)) {
        if($strand eq "+") {
            if($lastExon->{end}+3<=$substrateLen || $opt_b){$lastExon->{end}+=3}
            else{die "can't!  $lastExon->{end} $substrateLen"}}
	else
          {if($lastExon->{begin}-3>=0) {$lastExon->{begin}-=3}}
	undef $trans->{sequence};
	my $len=length $transSeq;
	undef $trans->getIthExon(0)->{sequence};
	if(!$opt_b) {$transSeq=$trans->loadTranscriptSeq($genomicSequence)}
	my $newLen=$opt_b ? 0 : length $transSeq;
	$stopLexeme=$opt_b ? "TAG" : substr($transSeq,length($transSeq)-3,3);
	if(!$isStop{$stopLexeme} && !$opt_f)
	  {
	    print "WARNING! $stopLexeme is not a stop codon; deleting transcript $transcriptId\n(exon type=$lastExonType)\n";
	    undef $transcripts->[$i];
	    next;
	  }
      }
    my $newGff=$trans->toGff();
    $gff.=$newGff;
    undef $trans->{sequence};
  }

open(OUT,">$outFilename") || die "can't write to file: $outFilename\n";
print OUT "$gff\n";
close(OUT);



