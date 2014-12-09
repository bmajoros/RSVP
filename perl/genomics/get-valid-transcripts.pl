#!/usr/bin/perl
use strict;
#use lib("/home/bmajoros/genomics/perl");
use FastaReader;
use FastaWriter;

my $usage="$0 <*.transcripts>   >  outfile";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my %stopCodons;
$stopCodons{"TAG"}=$stopCodons{"TGA"}=$stopCodons{"TAA"}=1;

my $reader=new FastaReader($filename);
my $writer=new FastaWriter;
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;

    my $startCodon=0;
    if($defline=~/startCodon=(\S+)/) {$startCodon=$1}
    $sequence=substr($sequence,$startCodon);

    # Ensure that the transcript starts with ATG
    next unless $sequence=~/^ATG/;

    # Ensure that the transcript ends with TAG/TGA/TAA
    my $len=length($sequence);
    my $numCodons=int($len/3);
    my $i;
    my $seq;
    for($i=0 ; $i<$numCodons ; ++$i)
      {
	my $codon=substr($sequence,$i*3,3);
	$seq.=$codon;
	last if $stopCodons{$codon};
      }
    next if $i>=$numCodons;

    # Output transcript
    $writer->addToFasta($defline,$seq,\*STDOUT);
  }


