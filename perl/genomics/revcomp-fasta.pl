#!/usr/bin/perl
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use FastaReader;
use FastaWriter;
use Translation;

#########################################################################
# Computes the reverse-complement of a fasta file and outputs the
# result as a new fasta file on STDOUT.
#
# bmajoros@tigr.org
#########################################################################

my $usage="$0 <*.fasta>\n";
die $usage unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new FastaReader($filename);
my $writer=new FastaWriter;
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    my $revcomp=Translation::reverseComplement(\$sequence);
    if($defline=~/^\s*>\s*(\S+)(.*)/)
      {
	my ($id,$rest)=($1,$2);
	$defline=">$id REVERSE-COMPLEMENTED $rest";
      }
    else
      {
	$defline=">REVERSE-COMPLEMENTED";
      }
    $writer->addToFasta($defline,$revcomp,\*STDOUT);
  }



