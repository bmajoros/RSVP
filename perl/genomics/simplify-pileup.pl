#!/usr/bin/perl
use strict;

my ($prevChr,$prevPos);
while(<STDIN>) {
  chomp;
  my @fields=split/\t/,$_;
  next unless @fields==6;
  my ($chr,$pos,$refBase,$count,$seq,$qual)=@fields;
  #my $DEBUG=$seq;
  $seq=~s/\$//g;
  $seq=~s/\^.//g;
  $seq=~s/[<>]+//g;
  while($seq=~/[+-](\d+)/) {
    my $len=$1;
    $seq=~s/[+-]$len.{$len}//;
  }
  my $newCount=length($seq);

  --$pos; ###
  print "$chr\t$pos\t$newCount\n";
}


