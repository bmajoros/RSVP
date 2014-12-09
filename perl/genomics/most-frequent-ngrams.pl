#!/usr/bin/perl
use strict;
use FastaReader;

$0=~/([^\/]+)\s*$/;
die "$1 <*.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
my @seqs;
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    push @seqs,$seq;
  }
my $nSeqs=@seqs;

for(my $n=2 ; $n<8 ; ++$n)
  {
    my %hash;
    for(my $i=0 ; $i<$nSeqs ; ++$i)
      {
	my $seq=$seqs[$i];
	my $len=length $seq;
	my $begin=0;
	my $end=$len-$n+1;
	for(my $j=$begin ; $j<$end ; ++$j)
	  {
	    my $ngram=substr($seq,$j,$n);
	    ++$hash{$ngram};
	  }
      }
    my ($largestCount,$bestNgram);
    my @keys=keys %hash;
    foreach my $key (@keys)
      {
	my $count=$hash{$key};
	if(!defined($largestCount) || $count>$largestCount)
	  {
	    $largestCount=$count;
	    $bestNgram=$key;
	  }
      }
    print "$largestCount $bestNgram\n";
  }



