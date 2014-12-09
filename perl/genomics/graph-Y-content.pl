#!/usr/bin/perl
use strict;
use FastaReader;

die "graph-Y-content.pl <*.fasta> <window-size>\n" unless @ARGV==2;
my ($fasta,$windowLen)=@ARGV;

my $reader=new FastaReader($fasta);
my (@CT,@N);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined($def);
    my $len=length($seq);
    for(my $i=0 ; $i<$len-$windowLen ; ++$i)
      {
	my $ct=0;
	for(my $j=$i ; $j<$i+$windowLen ; ++$j)
	  {
	    my $base=substr($seq,$j,1);
	    if($base=~/[CT]/) {++$ct}
	  }
	for(my $j=$i ; $j<$i+$windowLen ; ++$j)
	  {
	    $CT[$j]+=$ct;
	    $N[$j]+=$windowLen;
	  }
      }
  }
my $len=@CT;
for(my $i=0 ; $i<$len ; ++$i)
  {
    next unless $N[$i];
    my $p=$CT[$i]/$N[$i];
    print "$i $p\n";
  }
