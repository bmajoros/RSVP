#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;

$0=~/([^\/]+)$/;
my $usage="$1 <file-with-seq-IDs> <*.fasta>";
die "$usage\n" unless @ARGV==2;
my ($idFile,$fasta)=@ARGV;

my %keep;
open(IN,$idFile) || die "can't open file: $idFile\n";
while(<IN>)
  {
    if(/(\S+)/)
      {
	my $id=$1;
	$keep{$id}=1;
      }
  }
close(IN);

my $writer=new FastaWriter;
my $reader=new FastaReader($fasta);
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    $defline=~/^\s*>\s*(\S+)/;
    my $id=$1;
    if(!$keep{$id} && $id=~/(.*)[;\.@\[\(\/-]/) {$id=$1}
    if($keep{$id})
      {
	if($defline=~/\n/) {chop $defline}
	if($sequence=~/\n/) {chop $sequence}
	$writer->addToFasta($defline,$sequence,\*STDOUT);
      }
  }

