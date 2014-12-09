#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;
use FileSize;

my $name=ProgramName::get();
die "$name <1.fasta> <2.fasta>\n" unless @ARGV==2;
my ($file1,$file2)=@ARGV;

my $size1=FileSize::fileSize($file1);
my $size2=FileSize::fileSize($file2);

if($size1>$size2) {($file1,$file2)=($file2,$file1)}

my @seqs;
my $reader=new FastaReader($file1);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    $def=~/^>\s*(\S+)/;
    my $id=$1;
    $seq="\U$seq";
    push @seqs,[$seq,$id];
  }
my $n=@seqs;

my $numHits=0;
$reader=new FastaReader($file2);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless defined $def;
    $def=~/^>\s*(\S+)/;
    my $id=$1;
    $seq="\U$seq";
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $record=$seqs[$i];
	my ($otherSeq,$otherId)=@$record;
	if($seq eq $otherSeq)
	  {
	    print "$id\=$otherId\n";
	    ++$numHits;
	  }
      }
  }
print "$numHits total hits\n";

