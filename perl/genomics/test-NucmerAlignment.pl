#!/usr/bin/perl
use strict;
use NucmerAlignment;
use ProgramName;
use FastaReader;

my $name=ProgramName::get();
die "$name <delta-file> <1.fasta> <2.fasta>\n" unless @ARGV==3;
my ($deltaFile,$fasta1,$fasta2)=@ARGV;

my $alignment=new NucmerAlignment($deltaFile);
my $numHits=$alignment->getNumHits();

my $reader=new FastaReader($fasta1);
my ($def1,$seq1)=$reader->nextSequence();
$reader=new FastaReader($fasta2);
my ($def2,$seq2)=$reader->nextSequence();

for(my $i=0 ; $i<$numHits ; ++$i)
  {
    my $hit=$alignment->getIthHit($i);
    my $numIdent=$hit->{numIdent};
    my $cells=$hit->getCells();
    my $numCells=@$cells;
    my ($prevX,$prevY);
    my $numIdentities=0;
    my $alignLen=0;
    for(my $i=0 ; $i<$numCells ; ++$i)
      {
	my $cell=$cells->[$i];
	my ($x,$y)=@$cell;
	my $r1=substr($seq1,$x,1);
	my $r2=substr($seq2,$y,1);
	my $c=" ";
	if($x==$prevX) {$r1="|"}
	if($y==$prevY) {$r2="|"}
	if($r1 eq $r2) {$c="=";++$numIdentities}
	print "$r1 $c $r2          ($x,$y)\n";
	($prevX,$prevY)=($x,$y);
	++$alignLen;
      }
    my $pctIdent=int($numIdentities/$alignLen*100+5/9);
    print STDERR "$numIdentities identities = $pctIdent\% (expecting $numIdent)\n";
    exit;
  }



