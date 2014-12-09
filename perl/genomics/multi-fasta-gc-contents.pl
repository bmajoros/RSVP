#!/usr/bin/perl
use strict;
use FastaReader;

die "$0 <*.fasta>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $reader=new FastaReader($infile);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    $def=~/^\s*>\s*(\S+)/ || die "can't parse defline: $def\n";
    my $id=$1;
    my $gc=($seq=~s/([GC])/$1/g);
    my $gcat=($seq=~s/([GCAT])/$1/g);
    next unless $gcat>0;
    my $content=int(1000*$gc/$gcat+5/9)/10;
    print "$id\t$content\n";
  }
