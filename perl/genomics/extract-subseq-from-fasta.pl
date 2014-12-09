#!/usr/bin/perl
use strict;
use FastaReader;
use CompiledFasta;
use Getopt::Std;
our $opt_d;
getopts('d');

$0=~/([^\/]+)\s*$/;
my $usage="
$1 [-d] <*.fasta> <seq-id-or-dot(.)> <begin> <end>

   -d = preserve defline

   coordinates are zero-based/space-based
";
die "$usage\n" unless @ARGV==4;
my ($fasta,$seqId,$begin,$end)=@ARGV;
my $length=$end-$begin;

my $reader;
if($fasta=~/\.cof$/)
  {
    $reader=new CompiledFasta($fasta);
    my $subseq=$reader->load($begin,$end-$begin);
    print "$subseq\n";
  }
else
  {
    $reader=new FastaReader($fasta);
    while(1)
      {
	my ($defline,$seq)=$reader->nextSequence();
	last unless defined($defline);
	if($seqId ne ".")
	  {
	    $defline=~/$\s*>\s*([^\s;]+)/ || die;
	    my $id=$1;
	    next unless $id eq $seqId;
	  }
	my $subseq=substr($seq,$begin,$length);
	if(defined($opt_d)) {print "$defline"}
	print "$subseq\n";
	exit unless $seqId eq ".";
      }
  }

