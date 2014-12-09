#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffTranscriptReader;
use FastaReader;
use CodonIterator;

my $usage="$0 <*.gff> <*.fasta>";
die "$usage\n" unless @ARGV==2;
my ($filename,$fastaFilename)=@ARGV;

my %stopCodons=%{{"TAG"=>1,"TGA"=>1,"TAA"=>1}};

my $fastaReader=new FastaReader($fastaFilename);
my ($defline,$sequence)=$fastaReader->nextSequence();
my $substrateLength=length($sequence);

my $reader=new GffTranscriptReader();
my $transcriptArray=$reader->loadGFF($filename);
my $numTranscripts=@$transcriptArray;
my $substrate=($numTranscripts>0 ? $transcriptArray->[0]->{substrate} : ".");
print "$substrate\tvertex\tLEFT_TERMINUS\t0\t0\t.\t.\t.\n";
my $strand;

##########################################################################
# Rename exon types to things like "internal-exon", "initial-exon", etc...
##########################################################################
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    my $strand=$transcript->{strand};
    my $numExons=$transcript->numExons();
    last unless $transcript->getIthExon(0)->{type} eq "exon";
    $transcript->loadExonSequences(\$sequence);
    my $exon=$transcript->getIthExon(0);
    my $exonSeq=$exon->getSequence();
    my $firstThree=substr($exonSeq,0,3);
    if($firstThree eq "ATG") {$exon->{type}="initial-exon"}
    $exon=$transcript->getIthExon($numExons-1);
    $exonSeq=$exon->getSequence();
    my $lastThree=substr($exonSeq,length($exonSeq)-3,3);
    if($stopCodons{$lastThree})
      {
	if($numExons==1 && $exon->{type} eq "initial-exon")
	  {$exon->{type}="single-exon"}
	else
	  {$exon->{type}="final-exon"}
      }
    for(my $j=0 ; $j<$numExons ; ++$j)
      {
	my $exon=$transcript->getIthExon($j);
	if($exon->{type} eq "exon") {$exon->{type}="internal-exon"}
      }
  }

#######################################
# Adjust exon boundaries to exclude UTR
#######################################
goto SKIP;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    my $strand=$transcript->{strand};
    my $numExons=$transcript->numExons();
    my $startCodon=$transcript->{startCodon};
    if(!defined($startCodon)) {$startCodon=0} ### <---BAD IDEA!
    for(my $j=0 ; $j<$numExons ; ++$j)
      {
	my $exon=$transcript->getIthExon($j);
	my $length=$exon->getLength();
	if($length<=$startCodon)
	  {
	    $transcript->deleteExon($j);
	    --$numExons;
	    --$j;
	    $startCodon-=$length;
	  }
	else
	  {
	    if($strand eq "+")
	      {
		$exon->trimInitialPortion($startCodon);
		$transcript->{begin}=$exon->{begin};
	      }
	    else
	      {
		$exon->trimInitialPortion($startCodon);### ???
		$transcript->{end}=$exon->{end};
	      }
	    $exon->{type}=($numExons>1 ? "initial-exon" : "single-exon");
	    $transcript->{startCodon}=0;
	    last;
	  }
      }

    # Find in-frame stop codon
    my $codonIterator=new CodonIterator($transcript,\$sequence,\%stopCodons);
    my $stopCodonFound=0;
    while(my $codon=$codonIterator->nextCodon())
      {
	print "==>$codon->{triplet}\n";
	if($stopCodons{$codon->{triplet}})
	  {
	    my $exon=$codon->{exon};
	    my $coord=$codon->{absoluteCoord};
	    my $trimBases;
	    if($strand eq "+")
	      {$trimBases=$exon->{end}-$coord-3}
	    else
	      {$trimBases=$coord-$exon->{begin}-3}
	    $exon->trimFinalPortion($trimBases);
	    $exon->{type}=($exon->{order}==0 ? "single-exon" : "final-exon");
	    for(my $j=$numExons-1 ; $j>$exon->{order} ; --$j)
	      {$transcript->deleteExon($j)}
	    $stopCodonFound=1;
	    last;
	  }
	    #print STDERR "$codon->{triplet} is not a stop codon\n";
      }
    if(!$stopCodonFound) 
      {
	#die "stop codon not found for transcript $transcript->{transcriptId}"

	print STDERR "Warning!  No stop codon found for $transcript->{transcriptId} (skipping)\n";
	splice(@$transcriptArray,$i,1);
	--$i;
	--$numTranscripts;
      }

    $transcript->recomputeBoundaries();
  }
SKIP:

#####################################
# LEFTMOST INTERGENIC/INTRONIC REGION
#####################################
if($numTranscripts==0)
  {print "$substrate\tedge\tINTERGENIC\t1\t$substrateLength\t.\t.\t.\n"}
else
  {
    my $firstTranscript=$transcriptArray->[0];
    my $numExons=$firstTranscript->numExons();
    my $end=$firstTranscript->{begin};
    my $strand=$firstTranscript->getStrand();
    my $firstExon=($strand eq "+" ? $firstTranscript->getIthExon(0) :
		   $firstTranscript->getIthExon($numExons-1));
    if($firstExon->{begin}>0)
      {
	my $exonType=$firstExon->{type};
	if($strand eq "+" && 
	   ($exonType eq "internal-exon" || $exonType eq "final-exon") ||
	   $strand eq "-" &&
	   ($exonType eq "internal-exon" || $exonType eq "initial-exon"))
	  {#$end-=2;
	   print "$substrate\tedge\tINTRON\t1\t$end\t.\t$strand\t.\n"}
	else
	  {print "$substrate\tedge\tINTERGENIC\t1\t$end\t.\t.\t.\n"}
      }
  }

###########################################################
# HANDLE EACH TRANSCRIPT AND SURROUNDING INTERGENIC REGIONS
###########################################################
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    $strand=$transcript->{strand};
    my $numExons=$transcript->numExons();
    $substrate=$transcript->{substrate};

    # INTERGENIC REGION PRECEDING THIS TRANSCRIPT:
    if($i>0)
      {
	my $prevTranscript=$transcriptArray->[$i-1];
	my $begin=$prevTranscript->{end}+1;
	my $end=$transcript->{begin};
	print "$substrate\tedge\tINTERGENIC\t$begin\t$end\t.\t.\t.\n"
      }

    #------------------------------------------------------------------
    # PLUS STRAND
    #------------------------------------------------------------------
    if($strand eq "+")
      {
	# SINGLE EXON GENE
	if($numExons==1)
	  {
	    # ATG OR AG
	    my $exon=$transcript->getIthExon(0);
	    if($exon->{begin}>0)
	      {
		if($exon->{type} eq "final-exon" ||
		   $exon->{type} eq "internal-exon")
		  {printAG($exon)}
		else
		  {printATG($exon)}
	      }

	    # THE EXON
	    printExon($exon,$exon->{type});

	    # TAG OR GT
	    if($exon->{end}<$substrateLength)
	      {
		if($exon->{type} eq "initial-exon" ||
		   $exon->{type} eq "internal-exon")
		  {printGT($exon)}
		else
		  {printTAG($exon)}
	      }
	  }
	# MULTIPLE EXON GENE
	else
	  {
	    #####################
	    # HANDLE INITIAL EXON
	    #####################
	    my $initialExon=$transcript->getIthExon(0);
	    if($initialExon->{begin}>0)
	      {
		if($initialExon->{type} eq "initial-exon")
		  {printATG($initialExon)}
		else
		  {printAG($initialExon)}
	      }
	    printExon($initialExon,$initialExon->{type}); ### 12/7/02
	    printGT($initialExon); ### 12/7/02
	
	    #######################
	    # HANDLE INTERNAL EXONS
	    #######################
	    my $prevExon=$initialExon;
	    for(my $j=1 ; $j<$numExons-1 ; ++$j)
	      {
		my $exon=$transcript->getIthExon($j);
		printIntron($prevExon,$exon);
		printAG($exon);
		printExon($exon,$exon->{type});
		printGT($exon);
		$prevExon=$exon;
	      }
	
	    ###################
	    # HANDLE FINAL EXON
	    ###################
	    my $finalExon=$transcript->getIthExon($numExons-1);
	    printIntron($prevExon,$finalExon);
	    printAG($finalExon);
	    printExon($finalExon,$finalExon->{type});
	    if($finalExon->{end}<$substrateLength)
	      {
		if($finalExon->{type} eq "final-exon")
		  {printTAG($finalExon)}
		else
		  {printGT($finalExon)}
	      }
	  }
      }
    #------------------------------------------------------------------
    # MINUS STRAND
    #------------------------------------------------------------------
    else
      {
	# SINGLE EXON GENE
	if($numExons==1)
	  {
	    # TAG OR GT
	    my $exon=$transcript->getIthExon(0);
	    if($exon->{begin}>0)
	      {
		if($exon->{type} eq "single-exon" ||
		   $exon->{type} eq "final-exon")
		  {printTAG($exon)}
		else
		  {printGT($exon)}
	      }

	    # THE EXON
	    printExon($exon,$exon->{type});

	    # ATG OR AG
	    if($exon->{end}<$substrateLength)
	      {
		if($exon->{type} eq "single-exon" ||
		   $exon->{type} eq "initial-exon")
		  {printATG($exon)}
		else
		  {printAG($exon)}
	      }
	  }
	# MULTIPLE EXON GENE
	else
	  {
	    ###################
	    # HANDLE FIRST EXON
	    ###################
	    my $finalExon=$transcript->getIthExon($numExons-1);
	    if($finalExon->{begin}>0)
	      {
		if($finalExon->{type} eq "final-exon")
		  {printTAG($finalExon)}
		else
		  {printGT($finalExon)}
	      }
	    printExon($finalExon,$finalExon->{type}); ### 12/7/02
	    printAG($finalExon); ### 12/7/02

	    #######################
	    # HANDLE INTERNAL EXONS
	    #######################
	    my $prevExon=$finalExon;
	    for(my $j=$numExons-2 ; $j>0 ; --$j)
	      {
		my $exon=$transcript->getIthExon($j);
		printIntron($prevExon,$exon);
		printGT($exon);
		printExon($exon,$exon->{type});
		printAG($exon);
		$prevExon=$exon;
	      }
	
	    ##################
	    # HANDLE LAST EXON
	    ##################
	    my $initialExon=$transcript->getIthExon(0);
	    printIntron($prevExon,$initialExon);
	    printGT($initialExon);
	    printExon($initialExon,$initialExon->{type});
	    if($initialExon->{end}<$substrateLength)
	      {
		if($initialExon->{type} eq "initial-exon")
		  {printATG($initialExon)}
		else
		  {printAG($initialExon);}
	      }
	  }
      }
  }

######################################
# RIGHTMOST INTERGENIC/INTRONIC REGION
######################################
if($numTranscripts>0)
  {
    my $lastTranscript=$transcriptArray->[$numTranscripts-1];
    my $strand=$lastTranscript->getStrand();
    my $begin=$lastTranscript->{end}+1;
    my $numExons=$lastTranscript->numExons();
    my $lastExon=
      ($strand eq "+" ? $lastTranscript->getIthExon($numExons-1)
       : $lastTranscript->getIthExon(0));
    if($lastExon->{end}<$substrateLength)
      {
	my $exonType=$lastExon->{type};
	if($strand eq "+" && 
	   ($exonType eq "internal-exon" || $exonType eq "initial-exon") ||
	   $strand eq "-" && 
	   ($exonType eq "internal-exon" || $exonType eq "final-exon"))
	  {
	    #$begin+=2;
	    print "$substrate\tedge\tINTRON\t$begin\t$substrateLength".
	      "\t.\t$strand\t.\n";
	  }
	else
	  {print "$substrate\tedge\tINTERGENIC\t$begin\t$substrateLength".
	     "\t.\t.\t.\n"}
      }
  }
my $end=$substrateLength+1;
print "$substrate\tvertex\tRIGHT_TERMINUS\t$end\t$end\t.\t.\t.\n";


#-------------------------------------------------------------------
sub printATG
  {
    my ($initialExon)=@_;
    my ($atgBegin, $atgEnd);
    if($strand eq '+') {$atgBegin=$initialExon->{begin}+1}
    else {$atgBegin=$initialExon->{end}-2}
    $atgEnd=$atgBegin+2;
    print "$substrate\tvertex\tATG\t$atgBegin\t$atgEnd\t.\t$strand\t.\n";
  }

sub printTAG
  {
    my ($finalExon)=@_;
    my ($tagBegin, $tagEnd);
    if($strand eq "+") {$tagBegin=$finalExon->{end}-2}
    else {$tagBegin=$finalExon->{begin}+1}
    $tagEnd=$tagBegin+2;
    print "$substrate\tvertex\tTAG\t$tagBegin\t$tagEnd\t.\t$strand\t.\n";
  }
#-------------------------------------------------------------------
sub printGT
  {
    my ($exon)=@_;
    my ($gtBegin, $gtEnd);
    if($strand eq "+")
      {
	$gtBegin=$exon->{end}+1;
	$gtEnd=$gtBegin+1;
	print "$substrate\tvertex\tGT\t$gtBegin\t$gtEnd\t.\t$strand\t.\n";
      }
    else
      {
	$gtBegin=$exon->{begin}-1;
	$gtEnd=$gtBegin+1;
	print "$substrate\tvertex\tGT\t$gtBegin\t$gtEnd\t.\t$strand\t.\n";
      }
  }
#-------------------------------------------------------------------
sub printAG
  {
    my ($exon)=@_;
    my ($agBegin, $agEnd);
    if($strand eq "+")
      {
	$agBegin=$exon->{begin}-1;
	$agEnd=$agBegin+1;
	print "$substrate\tvertex\tAG\t$agBegin\t$agEnd\t.\t$strand\t.\n";
      }
    else
      {
	$agBegin=$exon->{end}+1;
	$agEnd=$agBegin+1;
	print "$substrate\tvertex\tAG\t$agBegin\t$agEnd\t.\t$strand\t.\n";
      }
  }
#-------------------------------------------------------------------
sub printExon
  {
    my ($exon,$exonType)=@_;
    my $begin=$exon->{begin}+1;
    my $end=$exon->{end};
    my $length=$end-$begin;
    my $frame=$exon->{frame};
    if($strand eq "+") {$frame=($frame+$length)%3}
    else {$frame=($frame-1)%3}
    my $transcript=$exon->getTranscript();
    my $transcriptId=$transcript->{transcriptId};
    $transcriptId=~s/;//g;
    my $source=$transcript->{source};
    my $extra="\ttransgrp=$transcriptId;\tsrc=$source;";
    print "$substrate\tedge\t$exonType\t$begin\t$end\t.\t$strand\t$frame\t$extra\n";
  }
#-------------------------------------------------------------------
sub printIntron
  {
    my ($prevExon,$thisExon)=@_;
    my $begin=$prevExon->{end}+1;#+3;
    my $end=$thisExon->{begin};#-2;
    print "$substrate\tedge\tINTRON\t$begin\t$end\t.\t$strand\t.\n";
  }
#-------------------------------------------------------------------

