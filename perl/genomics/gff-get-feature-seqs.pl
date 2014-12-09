#!/usr/bin/perl
use strict;
use GffReader;
use FastaReader;
use FastaWriter;
use CompiledFasta;
use Getopt::Std;
our $opt_c;
getopts('c:');

$0=~/([^\/]+)\s*$/;
die "\n$1 [-c id] <feature-label> <IN.gff> <IN.fasta> <OUT.fasta>
      where -c indicates that fasta file is 'compiled' and has given id
               (no defline nor whitespace characters)
\n"
  unless @ARGV==4;
my ($featureLabel,$inGff,$inFasta,$outFasta)=@ARGV;

my $gffReader=new GffReader;
my $features=$gffReader->loadGFF($inGff);
my %bySubstrate;
foreach my $feature (@$features)
  {
    my $substrate=$feature->getSubstrate();
    push @{$bySubstrate{$substrate}},$feature;
  }

open(OUT,">$outFasta") || die "Can't create file: $outFasta\n";
my $fastaWriter=new FastaWriter;
my $fastaReader=($opt_c ? new CompiledFasta($inFasta) :
		 new FastaReader($inFasta));
my $id=1;
my %substrateAlreadySeen;
while(1)
  {
    my ($def,$seq)=$opt_c ? (">$opt_c","") : $fastaReader->nextSequence();
    last unless $def;
    $def=~/>\s*(\S+)/;
    my $substrate=$1;
    next if $substrateAlreadySeen{$substrate};
    $substrateAlreadySeen{$substrate}=1;
    my $features=$bySubstrate{$substrate};
    next unless $features;
    foreach my $feature (@$features)
      {
	next unless $feature->getType() eq $featureLabel;
	my $begin=$feature->getBegin();
	my $end=$feature->getEnd();
	my $featureSeq=
	  $opt_c ?
	    $fastaReader->load($begin,$end-$begin) :
	    substr($seq,$begin,$end-$begin);
	my $defline=">$id $featureLabel on substrate $substrate";
	$fastaWriter->addToFasta($defline,$featureSeq,\*OUT);
	++$id;
	undef $featureSeq;
      }
    last if $opt_c;
  }
close(OUT);






