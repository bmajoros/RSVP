#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffTranscriptReader;
use FastaReader;
use Translation;

my $usage="$0 <genes.gff> <*.fasta>";
die "$usage\n" unless @ARGV==2;
my ($gffFile,$fastaFile)=@ARGV;

my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($gffFile);

my $fastaReader=new FastaReader($fastaFile);
my ($defline,$sequence)=$fastaReader->nextSequence();

my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $strand=$transcript->{strand};
    my $id=$transcript->{transcriptId};
    my $numExons=$transcript->numExons();
    for(my $j=0 ; $j<$numExons ; ++$j)
      {
	my $exon=$transcript->getIthExon($j);
	my ($begin,$end)=($exon->{begin},$exon->{end});
	if($j==0)
	  {# get start codon
	    if($strand eq "+")
	      {
		my $coord=$begin;
		my $startCodon=substr($sequence,$coord,3);
		++$coord;
		print "$coord START(+): $startCodon \t($id)\n";
	      }
	    else
	      {
		my $coord=$end-3;
		my $startCodon=substr($sequence,$coord,3);
		$startCodon=Translation::reverseComplement(\$startCodon);
		++$coord;
		print "$coord START(-): $startCodon \t($id)\n";
	      }
	  }
	else
	  {# get acceptor
	    if($strand eq "+")
	      {
		my $coord=$begin-2;
		my $acceptor=substr($sequence,$coord,2);
		++$coord;
		print "$coord ACCEPTOR(+): $acceptor \t($id)\n";
	      }
	    else
	      {
		my $coord=$end;
		my $acceptor=substr($sequence,$coord,2);
		$acceptor=Translation::reverseComplement(\$acceptor);
		++$coord;
		print "$coord ACCEPTOR(-): $acceptor \t($id)\n";
	      }
	  }

	if($j==$numExons-1)
	  {# get stop codon
	    if($strand eq "+")
	      {
		my $coord=$end-3;
		my $stopCodon=substr($sequence,$coord,3);
		++$coord;
		print "$coord STOP(+): $stopCodon \t($id)\n";
	      }
	    else
	      {
		my $coord=$begin;
		my $stopCodon=substr($sequence,$coord,3);
		$stopCodon=Translation::reverseComplement(\$stopCodon);
		++$coord;
		print "$coord STOP(-): $stopCodon \t($id)\n";
	      }
	  }
	else
	  {# get donor
	    if($strand eq "+")
	      {
		my $coord=$end;
		my $donor=substr($sequence,$coord,2);
		++$coord;
		print "$coord DONOR(+): $donor \t($id)\n";
	      }
	    else
	      {
		my $coord=$begin-2;
		my $donor=substr($sequence,$coord,2);
		$donor=Translation::reverseComplement(\$donor);
		++$coord;
		print "$coord DONOR(-): $donor \t($id)\n";
	      }
	  }
      }
  }



