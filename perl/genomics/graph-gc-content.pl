#!/usr/bin/perl
use strict;
use FastaReader;
$|=1;

die "graph-gc-content.pl <*.fasta> <window-size>\n" unless @ARGV==2;
my ($fasta,$windowLen)=@ARGV;

my $INCREMENT=$windowLen/100;

my $reader=new FastaReader($fasta);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined($def);
    my $len=length($seq);
    my (@GC,@N);
    for(my $i=0 ; $i<$len-$windowLen ; $i+=$INCREMENT)
      {
	my $gc=0;
	#for(my $j=$i ; $j<$i+$windowLen ; ++$j)
	#  {
	#    my $base=substr($seq,$j,1);
	#    if($base=~/[GC]/) {++$gc}
	#  }
	my $wnd=substr($seq,$i,$windowLen);
	$gc=($wnd=~s/([GC])/$1/g);
	$GC[$i]+=$gc;
	$N[$i]+=$windowLen;
      }
    for(my $i=0 ; $i<$len ; $i+=$INCREMENT)
      {
	next unless $N[$i];
	my $p=$GC[$i]/$N[$i];
	print "$i $p\n";
      }
last;
  }
