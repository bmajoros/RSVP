#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use SummaryStats;
use Translation;
use FastaWriter;
use FastaReader;
use FileHandle;
use GffReader;
use Transcript;
use CodonIterator;
use Codon;
$|=1;

my $usage="$0 <*.gff> <*.seq>";
die "$usage\n" unless @ARGV==2;
my ($gffFilename,$seqFilename)=@ARGV;

# Read the sequence from the FASTA file
my $fastaReader=new FastaReader($seqFilename);
my ($defline,$axisSequence)=$fastaReader->nextSequence();

# Read the feature coordinates from the GFF file
my $gffReader=new GffReader;
my $transcripts=$gffReader->loadGFF($gffFilename);

my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $trans=$transcripts->[$i];
    my $transcriptId=$trans->{transcriptId};
    my $strand=$trans->{strand};
    print "transcript $transcriptId: ($strand)\n";

    my $iterator=new CodonIterator($trans,\$axisSequence);
    my $codons=$iterator->getAllCodons();
    my $numCodons=@$codons;
    my $firstCodon=$codons->[0];
    my $lastCodon=$codons->[$numCodons-1];

    my $triplet=$firstCodon->{triplet};
    my $exonOrder=$firstCodon->{exon}->{order};
    my $absolute=$firstCodon->{absolute};
    my $relative=$firstCodon->{relative};
    my $isInterrupted=$firstCodon->{isInterrupted} || 0;
    my $basesInExon=$firstCodon->{basesInExon};
    my $exonLen=$firstCodon->{exon}->getLength();
    print "\tSTART: $triplet I=$isInterrupted B=$basesInExon O=$exonOrder\tA=$absolute\tR=$relative\tEL=$exonLen\n";
    
    $triplet=$lastCodon->{triplet};
    $exonOrder=$lastCodon->{exon}->{order};
    $absolute=$lastCodon->{absolute};
    $relative=$lastCodon->{relative};
    $isInterrupted=$lastCodon->{isInterrupted} || 0;
    $basesInExon=$lastCodon->{basesInExon};
    $exonLen=$lastCodon->{exon}->getLength();
    print "\tEND: $triplet I=$isInterrupted B=$basesInExon O=$exonOrder\tA=$absolute\tR=$relative\tEL=$exonLen\n";
    

    next;


    while(1)
      {
	my $codon=$iterator->nextCodon();
	last unless defined $codon;
	my $triplet=$codon->{triplet};
	my $exonOrder=$codon->{exon}->{order};
	die "no exon" unless defined $codon->{exon};
	if(!defined $codon->{exon}->{order})
	  {
	    my $exon=$codon->{exon};
	    my @keys=keys %$exon;
	    print "Exon $exon:\n";
	    foreach my $key (@keys)
	      {
		my $value=$exon->{$key};
		print "\t$key -> $value\n";
	      }
	    die "no order for exon";
	  }
	my $absolute=$codon->{absolute};
	my $relative=$codon->{relative};
	my $isInterrupted=$codon->{isInterrupted} || 0;
	my $basesInExon=$codon->{basesInExon};
	my $exonLen=$codon->{exon}->{end}-$codon->{exon}->{begin};###debugging
	print "\t$triplet I=$isInterrupted B=$basesInExon O=$exonOrder\tA=$absolute\tR=$relative\tEL=$exonLen\n";
      }

  }
