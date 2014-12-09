#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <infile> <outfile>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't create $outfile";
my $writer=new FastaWriter;

my @seqs;
my $reader=new FastaReader($infile);
while(1)
  {
    my ($defline,$seq)=$reader->nextSequence();
    last unless defined $defline;
    push @seqs,[\$defline,\$seq];
  }

my $n=@seqs;
for(my $i=0 ; $i<$n ; ++$i)
{
    my $j=int(rand($n));
    my $a=$seqs[$i];
    my $b=$seqs[$j];
    $seqs[$i]=$b;
    $seqs[$j]=$a;
}

for(my $i=0 ; $i<$n ; ++$i)
{
    my $a=$seqs[$i];
    my ($def,$seq)=@$a;
    $writer->addToFasta($$def,$$seq,\*OUT);
}
close(OUT);



