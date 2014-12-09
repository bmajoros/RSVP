#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.fasta> TAG,TGA,TAA [transcript-id]";
die "$usage\n" unless @ARGV==2 || @ARGV==3;
my ($filename,$stops,$transcriptId)=@ARGV;

my @stops=split/,/,$stops;
my %stopCodons;#=%{{TAG=>1,TGA=>1,TAA=>1}};
foreach my $stop (@stops) {$stopCodons{$stop}=1}

my $fastaReader=new FastaReader($filename);
my $phase=0;
while(1)
  {
    my ($defline,$sequence)=$fastaReader->nextSequence();
    if(!defined($defline)) {last}
    if(length($transcriptId) && $defline=~/^\s*>(\S+)/)
      {
	my $id=$1;
	if($id=~/(.*);/) {$id=$1}
	next unless($transcriptId eq $id);
      }
    if($defline=~/frame=(\d)/) {$phase=$1}
    print "$defline";
    my $discardBases=(3-$phase)%3;
    my $length=length($sequence);
    $sequence=substr($sequence,$discardBases,$length-$discardBases);
    $length-=$discardBases;
    if($defline=~/startCodon=(\d+)/)
      {$sequence=substr($sequence,$1,$length)}
    my $numCodons=int($length/3);
    for(my $i=0 ; $i<$numCodons ; ++$i)
      {
	my $codon=substr($sequence,$i*3,3);
	print "$codon ";
#	if($stopCodons{$codon}) {last}
	if($stopCodons{$codon} && $i<$numCodons-1) 
	  {die "in-frame stop codon $codon at pos=$i\n"}
	if(($i+1)%15==0) {print "\n"}
      }
    print "\n";
  }




