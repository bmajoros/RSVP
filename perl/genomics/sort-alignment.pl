#!/usr/bin/perl
use strict;
use ProgramName;
use FastaReader;
use FastaWriter;

my $name=ProgramName::get();
die "$name <in.fasta> <out.fasta> <tree.newick>\n" unless @ARGV==3;
my ($infile,$outfile,$treeFile)=@ARGV;

my @species;
open(IN,$treeFile) || die "can't open $treeFile\n";
my $tree;
while(<IN>) {
  if(/tree\s+\S+\s*=\s*(\S.*\S)/) { $tree=$1 }
  elsif(/\(\S.*\S\)/) { $tree=$1 }
}
$tree=~s/\s+//g;
$tree=~s/:[\d\.]+//g;
my @fields=split/[\(\),;]+/,$tree;
foreach my $field (@fields) {
  if($field=~/\S/) { push @species,$field }
}
close(IN);

my %seqs;
my $reader=new FastaReader($infile);
while(1) {
  my ($def,$seq)=$reader->nextSequence();
  last unless length($seq)>0;
  $def=~/^>(\S+)/ || die;
  my $ID=$1;
  $seqs{$ID}.=$seq;
}

my $writer=new FastaWriter;
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
foreach my $ID (@species) {
  my $seq=$seqs{$ID};
  next unless $seq;
  $seq="\U$seq";
  $seq=~s/\s+//g;
  $seq=~s/[^ACGT]/-/g;
  $writer->addToFasta(">$ID",$seq,\*OUT);
}
close(OUT);
