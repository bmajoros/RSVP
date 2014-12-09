package Transcript;
use strict;
use Carp;
use Exon;
use Translation;
use CodonIterator;
use CompiledFasta;
use Carp;

######################################################################
#
# Transcript.pm
#
# bmajoros@tigr.org 7/8/2002
#
#
# Attributes:
#   substrate : scaffold or chromosome name
#   source : name of entity that predicted or curated this transcript
#   startCodon : index into (spliced) transcript sequence of ATG,
#                regardless of strand
#   startCodonAbsolute : absolute coordinates of start codon, 
#                        relative to genomic axis
#   strand : + or -
#   exons : pointer to array of Exons
#   begin : begin coordinate of leftmost exon (zero based)
#   end : end coordinate of rightmost exon (one past end)
#   sequence : NOT LOADED BY DEFAULT!
#   transcriptId : identifier
#   stopCodons : hash table of stop codons (strings)
#   geneId : identifier of gene to which this transcript belongs
#   gene : a Gene object to which this transcript belongs
# Methods:
#   $transcript=new Transcript($id,$strand);
#   $transcript->addExon($exon);
#   $copy=$transcript->copy();
#   $bool=$transcript->areExonTypesSet();
#   $transcript->setExonTypes();
#   $success=$transcript->loadExonSequences(\$axisSequence);
#   $seq=$transcript->loadTranscriptSeq(\$axisSequence);
#   $success=$transcript->loadExonSequencesCF($compiledFasta);
#   $seq=$transcript->loadTranscriptSeqCF($compiledFasta);
#   $bool=$transcript->equals($other);
#   $bool=$transcript->overlaps($otherTranscript);
#   $bool=$transcript->overlaps($begin,$end);
#   $bool=$transcript->isPartial();
#   $bool=$transcript->isContainedWithin($begin,$end);
#   $bool=$transcript->exonOverlapsExon($exon);
#   $len=$transcript->getLength(); # sum of exon sizes
#   $len=$transcript->getExtent(); # end-begin
#   $n=$transcript->numExons();
#   $exon=$transcript->getIthExon($i);
#   $transcript->deleteExon($index);
#   $transcript->deleteExonRef($exon);
#   $transcript->recomputeBoundaries();# for after trimming 1st & last exons
#   $transcript->getSubstrate();
#   $transcript->getSource();
#   $gff=$transcript->toGff();
#   $id=$transcript->getID();
#   $id=$transcript->getTranscriptId();
#   $id=$transcript->getGeneId();
#   $transcript->setGeneId($id);
#   $transcript->setTranscriptId($id);
#   $transcript->setSubstrate($substrate);
#   $transcript->setSource($source);
#   $begin=$transcript->getBegin();
#   $end=$transcript->getEnd();
#   $transcript->setBegin($x);
#   $transcript->setEnd($x);
#   $strand=$transcript->getStrand();
#   if($transcript->isWellFormed(\$sequence)) ...   # See notes in the sub
#   $transcript->trimUTR();
#   $transcript->getScore();
#   my $introns=$transcript->getIntrons();
#   my $nextExon=$transcript->getSuccessor($thisExon);
#   $transcript->shiftCoords($delta);
#   $transcript->reverseComplement($seqLen);
#   $transcript->setStopCodons({TGA=>1,TAA=>1,TAG=>1});
#   $g=$transcript->getGene();
#   $transcript->setGene($g);
#   $genomicCoord=$transcript->mapToGenome($transcriptCoord);
#   $transcriptCoord=$transcript->mapToTranscript($genomicCoord);
#   $exon=$transcript->getExonContaining($genomicCoord);
#   $transcript->reverseExonOrder();
# Private methods:
#   $transcript->adjustOrders();
#   $transcript->sortExons();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class,$id,$strand)=@_;

  my $self=
    {
     transcriptId=>$id,
     strand=>$strand,
     exons=>[],
     stopCodons=>{TAG=>1,TGA=>1,TAA=>1},
    };
  bless $self,$class;

  return $self;
}
#------------------------------------------------------
#   $bool=$transcript->equals($other);
sub equals
  {
    my ($this,$other)=@_;
    if($this->getSource() ne $other->getSource()) {return 0}
    if($this->getStrand() ne $other->getStrand()) {return 0}
    if($this->getBegin()  != $other->getBegin())  {return 0}
    if($this->getEnd()    != $other->getEnd())    {return 0}
    my $n=$this->numExons();
    my $m=$other->numExons();
    if($n!=$m) {return 0}
    for(my $i=0 ; $i<$m ; ++$i)
      {
	my $thisExon=$this->getIthExon($i);
	my $thatExon=$other->getIthExon($i);
	if($thisExon->getBegin() != $thatExon->getBegin()) {return 0}
	if($thisExon->getEnd()   != $thatExon->getEnd())   {return 0}
      }
    return 1;
  }
#------------------------------------------------------
sub compStrand
  {
    my ($strand)=@_;
    if($strand eq "+") {return "-"}
    if($strand eq "-") {return "+"}
    if($strand eq ".") {return "."}
    die "Unknown strand \"$strand\" in $0";
  }
#---------------------------------------------------------------------
#   $transcript->setTranscriptId($id);
sub setTranscriptId
  {
    my ($transcript,$newId)=@_;
    $transcript->{transcriptId}=$newId;
  }
#---------------------------------------------------------------------
#   $transcript->reverseComplement($seqLen);
sub reverseComplement
  {
    my ($self,$seqLen)=@_;
    my $begin=$self->getBegin();
    my $end=$self->getEnd();
    $self->{begin}=$seqLen-$end;
    $self->{end}=$seqLen-$begin;
    $self->{strand}=compStrand($self->{strand});
    my $exons=$self->{exons};
    foreach my $exon (@$exons)
      {
	$exon->reverseComplement($seqLen);
      }
  }
#---------------------------------------------------------------------
#   $bool=$transcript->exonOverlapsExon($exon);
sub exonOverlapsExon
  {
    my ($self,$exon)=@_;
    my $exons=$self->{exons};
    my $n=@$exons;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $myExon=$exons->[$i];
	if($exon->overlaps($myExon)) {return 1}
      }
    return 0;
  }
#---------------------------------------------------------------------
#   $exon->shiftCoords($delta);
sub shiftCoords
  {
    my ($self,$delta)=@_;
    my $exons=$self->{exons};
    foreach my $exon (@$exons)
      {
	$exon->shiftCoords($delta);
      }
    $self->{begin}+=$delta;
    $self->{end}+=$delta;
  }
#---------------------------------------------------------------------
#   $success=$transcript->loadExonSequences($axisSequenceRef);
sub loadExonSequences
  {
    my ($self,$axisSequenceRef)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $strand=$self->{strand};
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my $start=$exon->{begin};
	my $length=$exon->{end}-$start;
	my $exonSeq=substr($$axisSequenceRef,$start,$length);
	if(length($exonSeq)!=$length)
	  {
	    confess "start=$start, length=$length, but substrate $self->{substrate} ends at ".
	      length($$axisSequenceRef);
	   }
	if($strand eq "-")
	  {$exonSeq=Translation::reverseComplement(\$exonSeq)}
	$exon->{sequence}=$exonSeq;
      }
    return 1;
  }
#---------------------------------------------------------------------
#   $seq=$transcript->loadTranscriptSeq($axisSequenceRef);
sub loadTranscriptSeq
  {
    my ($self,$axisSequenceRef)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $firstExon=$exons->[0];
    $self->loadExonSequences($axisSequenceRef)
      unless defined $firstExon->{sequence};
    my $sequence;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	$sequence.=$exon->{sequence};
      }
    $self->{sequence}=$sequence;
    return $sequence;
  }
#---------------------------------------------------------------------
#   $success=$transcript->loadExonSequencesCF($compiledFasta);
sub loadExonSequencesCF
  {
    my ($self,$compiledFasta)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $strand=$self->{strand};
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my $start=$exon->{begin};
	my $length=$exon->{end}-$start;
	my $exonSeq=$compiledFasta->load($start,$length);
	$exonSeq="\U$exonSeq";
	if(length($exonSeq)!=$length)
	  {
	    confess "Error loading exon sequence in Transcript::loadExonSequencesCF()";
	   }
	if($strand eq "-")
	  {$exonSeq=Translation::reverseComplement(\$exonSeq)}
	$exon->{sequence}=$exonSeq;
      }
    return 1;
  }
#---------------------------------------------------------------------
#   $seq=$transcript->loadTranscriptSeqCF($compiledFasta);
sub loadTranscriptSeqCF
  {
    my ($self,$compiledFasta)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $firstExon=$exons->[0];
    $self->loadExonSequencesCF($compiledFasta)
      unless defined $firstExon->{sequence};
    my $sequence;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	$sequence.=$exon->{sequence};
      }
    $self->{sequence}=$sequence;
    return $sequence;
  }
#---------------------------------------------------------------------
#   $bool=$transcript->isContainedWithin($begin,$end);
sub isContainedWithin
  {
    my($this,$begin,$end)=@_;
    return $this->{begin}>=$begin && $this->{end}<=$end;
  }
#---------------------------------------------------------------------
#   $bool=$transcript->overlaps($begin,$end);
sub overlaps
  {
    my($this,$begin,$end)=@_;
    return $this->{begin}<$end && $begin<$this->{end};
  }
#---------------------------------------------------------------------
#   $bool=$transcript->overlaps($otherTranscript)
sub overlaps
  {
    my($this,$otherTranscript)=@_;
    return $this->{begin}<$otherTranscript->{end} &&
      $otherTranscript->{begin}<$this->{end};
  }
#---------------------------------------------------------------------
#   $len=$transcript->getLength(); # sum of exon sizes
sub getLength
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $nExons=@$exons;
    my $length=0;
    for(my $i=0 ; $i<$nExons ; ++$i)
      {
	my $exon=$exons->[$i];
	$length+=$exon->getLength();
      }
    return $length;
  }
#---------------------------------------------------------------------
#   $len=$transcript->getExtent(); # end-begin
sub getExtent
  {
    my ($self)=@_;
    return $self->{end}-$self->{begin};
  }
#---------------------------------------------------------------------
#   $n=$transcripts->numExons();
sub numExons
  {
    my ($self)=@_;
    return 0+@{$self->{exons}};
  }
#---------------------------------------------------------------------
#   $exon=$transcript->getIthExon($i);
sub getIthExon
  {
    my ($self,$i)=@_;
    return $self->{exons}->[$i];
  }
#---------------------------------------------------------------------
#   $transcript->deleteExon($i);
sub deleteExon
  {
    my ($self,$index)=@_;
    my $exons=$self->{exons};
    splice(@$exons,$index,1);
    $self->adjustOrders();
  }
#---------------------------------------------------------------------
#   $transcript->sortExons();
sub sortExons
  {
    my ($self)=@_;
    if($self->{strand} eq "+")
      {@{$self->{exons}}=sort {$a->{begin}<=>$b->{begin}} @{$self->{exons}}}
    else
      {@{$self->{exons}}=sort {$b->{begin}<=>$a->{begin}} @{$self->{exons}}}
  }
#---------------------------------------------------------------------
#   $transcript->adjustOrders();
sub adjustOrders
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
#    confess($numExons) unless $numExons>0;
    return unless $numExons>0;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {$exons->[$i]->{order}=$i}
    if($self->{strand} eq "+")
      {
	$self->{begin}=$exons->[0]->{begin};
	$self->{end}=$exons->[$numExons-1]->{end};
      }
    else
      {
	$self->{begin}=$exons->[$numExons-1]->{begin};
	$self->{end}=$exons->[0]->{end};
      }
  }
#---------------------------------------------------------------------
#   $bool=$transcript->areExonTypesSet();
sub areExonTypesSet
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my %validExonTypes=
      %{{"single-exon"=>1,
	 "initial-exon"=>1,
         "internal-exon"=>1,
         "final-exon"=>1}};
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $type=$exons->[$i]->getType();
	if(!$validExonTypes{$type}) {return 0}
      }
    return 1;
  }
#---------------------------------------------------------------------
#   $transcript->setExonTypes();
sub setExonTypes
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my %validExonTypes=
      %{{"single-exon"=>1,
	 "initial-exon"=>1,
         "internal-exon"=>1,
         "final-exon"=>1}};
    if($numExons==1)
      {$exons->[0]->setType("single-exon")
	 unless $validExonTypes{$exons->[0]->getType()}}
    else
      {
	for(my $i=1 ; $i<$numExons-1 ; ++$i)
	  {
	    $exons->[$i]->setType("internal-exon")
	      unless $validExonTypes{$exons->[$i]->getType()}
	  }
	$exons->[0]->setType("initial-exon");
	$exons->[$numExons-1]->setType("final-exon");
      }
  }
#---------------------------------------------------------------------
#   $transcript->deleteExonRef($exon);
sub deleteExonRef
  {
    my ($self,$victim)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $i;
    for($i=0 ; $i<$numExons ; ++$i)
      {
	my $thisExon=$exons->[$i];
	if($thisExon==$victim) {last}
      }
    if($i>=$numExons)
      {die "Can't find exon $victim in Transcript::deleteExon()"}
    $self->deleteExon($i);
  }
#---------------------------------------------------------------------
#   $transcript->recomputeBoundaries(); # for after trimming 1st & last exons
sub recomputeBoundaries
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $firstExon=$exons->[0];
    my $lastExon=$exons->[$numExons-1];
    my $strand=$self->{strand};
    if($strand eq "+")
      {
	$self->{begin}=$firstExon->{begin};
	$self->{end}=$lastExon->{end};
      }
    else
      {
	$self->{begin}=$lastExon->{begin};
	$self->{end}=$firstExon->{end};
      }
  }
#---------------------------------------------------------------------
#   $transcript->addExon($exon);
sub addExon
  {
    my ($self,$exon)=@_;
    my $strand=$exon->getStrand();
    my $exons=$self->{exons};
    push @$exons,$exon;
    if($strand eq "+")
      {@$exons=sort {$a->{begin} <=> $b->{begin}} @$exons}
    else
      {@$exons=sort {$b->{begin} <=> $a->{begin}} @$exons}
    $self->adjustOrders();
  }
#---------------------------------------------------------------------
#   $transcript->getSubstrate();
sub getSubstrate
  {
    my ($self)=@_;
    return $self->{substrate};
  }
#---------------------------------------------------------------------
#   $transcript->getSource();
sub getSource
  {
    my ($self)=@_;
    return $self->{source};
  }
#---------------------------------------------------------------------
#   my $gff=$transcript->toGff();
sub toGff
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    my $gff;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my $exonGff=$exon->toGff();
	$gff.=$exonGff;
      }
    return $gff;
  }
#---------------------------------------------------------------------
#   $id=$transcript->getTranscriptId();
sub getTranscriptId
  {
    my ($self)=@_;
    return $self->{transcriptId};
  }
#---------------------------------------------------------------------
#   $id=$transcript->getGeneId();
sub getGeneId
  {
    my ($self)=@_;
    return $self->{geneId};
  }
#---------------------------------------------------------------------
#   $transcript->setGeneId($id);
sub setGeneId
  {
    my ($self,$id)=@_;
    $self->{geneId}=$id;
  }
#---------------------------------------------------------------------
#   $id=$transcript->getID();
sub getID
  {
    my ($self)=@_;
    return $self->{transcriptId};
  }
#---------------------------------------------------------------------
#   $transcript->setSubstrate($substrate);
sub setSubstrate
  {
    my ($self,$substrate)=@_;
    $self->{substrate}=$substrate;
  }
#---------------------------------------------------------------------
#   $transcript->setSource($source);
sub setSource
  {
    my ($self,$source)=@_;
    $self->{source}=$source;
  }
#---------------------------------------------------------------------
#   $begin=$transcript->getBegin();
sub getBegin
  {
    my ($self)=@_;
    return $self->{begin};
  }
#---------------------------------------------------------------------
#   $end=$transcript->getEnd();
sub getEnd
  {
    my ($self)=@_;
    return $self->{end};
  }
#---------------------------------------------------------------------
#   $strand=$transcript->getStrand();
sub getStrand
  {
    my ($self)=@_;
    return $self->{strand};
  }
#---------------------------------------------------------------------
#   if($transcript->isWellFormed(\$sequence)) ...
#
#   This procedure iterates through the codons of this transcript,
#   starting at the start codon (attribute startCodon specifies this
#   offset within the transcript, not counting intron bases), and
#   continuing until either an in-frame stop codon is encountered,
#   or the end of the transcript is reached.  The transcript is
#   considered well-formed only if a stop-codon is encountered.
#
sub isWellFormed
  {
    my ($self,$seq)=@_;
    my $stopCodons=$self->{stopCodons};

    # 1. Check whether any exons overlap each other
    my $exons=$self->{exons};
    my $numExons=@$exons;
    for(my $i=1 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my $prevExon=$exons->[$i-1];
	if($exon->overlaps($prevExon)) {return 0}
      }

    # 2. Check that there is an in-frame stop-codon
    my $iterator=new CodonIterator($self,$seq,$stopCodons);
    my $codons=$iterator->getAllCodons();
    my $n=@$codons;
    return 0 unless $n>0;
    my $lastCodon=$codons->[$n-1];
    return $stopCodons->{$lastCodon->{triplet}};
  }
#---------------------------------------------------------------------
#   $transcript->trimUTR(\$axisSequence);
sub trimUTR
  {
    my ($self,$axisSequenceRef)=@_;
    $self->adjustOrders();
    my $stopCodons=$self->{stopCodons};
    my $sequence=$self->{sequence};

    my $strand=$self->{strand};
    my $numExons=$self->numExons();
    my $startCodon=$self->{startCodon};

    if(!defined($startCodon)) 
      {die "can't trim UTR, because startCodon is not set"}
    for(my $j=0 ; $j<$numExons ; ++$j)
      {
	my $exon=$self->getIthExon($j);
	my $length=$exon->getLength();
	if($length<=$startCodon)
	  {
	    $self->deleteExon($j);
	    --$numExons;
	    --$j;
	    $startCodon-=$length;
	    $self->adjustOrders(); ### 4/1/03
	  }
	else
	  {
	    if($strand eq "+")
	      {
		$exon->trimInitialPortion($startCodon);
		$self->{begin}=$exon->{begin};
	      }
	    else
	      {
		$exon->trimInitialPortion($startCodon);### ???
		$self->{end}=$exon->{end};
	      }
	    $exon->{type}=($numExons>1 ? "initial-exon" : "single-exon");
	    $self->{startCodon}=0;
	    last;
	  }
      }

    # Find in-frame stop codon
    my $codonIterator=
      new CodonIterator($self,$axisSequenceRef,$stopCodons);
    my $stopCodonFound=0;
    while(my $codon=$codonIterator->nextCodon())
      {
	if($stopCodons->{$codon->{triplet}})
	  {
	    my $exon=$codon->{exon};
	    my $coord=$codon->{absoluteCoord};
	    my $trimBases;
	    if($strand eq "+")
	      {$trimBases=$exon->{end}-$coord-3}
	    else
	      {$trimBases=$coord-$exon->{begin}-3}

	    if($trimBases>=0)
	      {
		$exon->trimFinalPortion($trimBases);
		$exon->{type}=
		  ($exon->{order}==0 ? "single-exon" : "final-exon");
		for(my $j=$numExons-1 ; $j>$exon->{order} ; --$j)
		  {$self->deleteExon($j)}
		$stopCodonFound=1;
		last;
	      }
	    else # codon is interrupted; trim the next exon
	      {
		my $nextExon=$self->getSuccessor($exon);
		if(!defined($nextExon)) {die "exon successor not found"}
		$nextExon->trimFinalPortion($nextExon->getLength()+$trimBases);
		$nextExon->{type}="final-exon";
		for(my $j=$numExons-1 ; $j>$nextExon->{order} ; --$j)
		  {$self->deleteExon($j)}
		$stopCodonFound=1;
		last;		
	      }
	  }
      }
    if(!$stopCodonFound)
      {
	### sometimes the GFF coords don't include the stop codon...
	my $numExons=$self->numExons();
	my $lastExon=$self->getIthExon($numExons-1);
	my $lastExonEnd=$lastExon->getEnd();
	my $seq=$axisSequenceRef;
	if($strand eq "+")
	  {
	    my $stopCodonBegin=$lastExonEnd;
	    my $stopCodon=substr($$seq,$stopCodonBegin,3);
	    if($stopCodon ne "TAG" && $stopCodon ne "TAA" 
	       && $stopCodon ne "TGA")
	      {
		print "Warning!  No stop codon found for $self->{transcriptId}, $self->{strand} strand, unable to trim UTR\n";
	      }
	  }
	else # $strand eq "-"
	  {
	    my $stopCodonBegin=$lastExon->getBegin()-3;
	    my $stopCodon=substr($$seq,$stopCodonBegin,3);
	    $stopCodon=Translation::reverseComplement(\$stopCodon);
	    if($stopCodon ne "TAG" && $stopCodon ne "TAA" 
	       && $stopCodon ne "TGA")
	      {
		print "Warning!  No stop codon found for $self->{transcriptId}, $self->{strand} strand, unable to trim UTR\n";
	      }
	  }
	#print "Warning!  No stop codon found for $self->{transcriptId} (skipping), $self->{strand} strand\n";
      }
    $self->recomputeBoundaries();
  }
#---------------------------------------------------------------------
#   $transcript->getScore();
sub getScore
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $n=@$exons;
    my $score=0;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $exon=$exons->[$i];
	my $exonScore=$exon->getScore();
	$score+=$exonScore unless $exonScore eq ".";
      }
    return $score;
  }
#---------------------------------------------------------------------
#   my $introns=$transcript->getIntrons();
sub getIntrons
  {
    my ($self)=@_;
    my $numExons=$self->numExons();
    my $strand=$self->{strand};
    my @introns;
    my ($lastExonEnd,$lastExonFrame);
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$self->getIthExon($i);
	if(defined($lastExonEnd))
	  {
	    if($strand eq "+")
	      {push @introns,[$lastExonEnd,$exon->getBegin(),$exon->getFrame()]}
	    else
	      {push @introns,[$exon->getEnd(),$lastExonEnd,$lastExonFrame]}
	  }
	$lastExonEnd=($strand eq "+" ? $exon->getEnd() : $exon->getBegin());
	$lastExonFrame=$exon->getFrame();
      }
    return \@introns;
  }
#---------------------------------------------------------------------
#   my $nextExon=$transcript->getSuccessor($thisExon);
sub getSuccessor
  {
    my ($self,$targetExon)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    for(my $i=0 ; $i<$numExons-1 ; ++$i)
      {
	my $exon=$exons->[$i];
	if($exon==$targetExon) {return $exons->[$i+1]}
      }
    return undef;
  }
#---------------------------------------------------------------------
#   $transcript->setStopCodons({TGA=>1,TAA=>1,TAG=>1});
sub setStopCodons
  {
    my ($self,$stopCodons)=@_;
    $self->{stopCodons}=$stopCodons;
  }
#---------------------------------------------------------------------
#   $bool=$transcript->isPartial();
sub isPartial
  {
    my ($self)=@_;
    my $exons=$self->{exons};
    my $numExons=@$exons;
    if($numExons==1) {return $exons->[0]->getType() ne "single-exon"}
    return 
      $exons->[0]->getType() ne "initial-exon" ||
      $exons->[$numExons-1]->getType() ne "final-exon";
  }
#---------------------------------------------------------------------
#   $g=$transcript->getGene();
sub getGene
  {
    my ($self)=@_;
    return $self->{gene};
  }
#---------------------------------------------------------------------
#   $transcript->setGene($g);
sub setGene
  {
    my ($self,$g)=@_;
    $self->{gene}=$g;
  }
#---------------------------------------------------------------------
#   $transcript->setBegin($x);
sub setBegin
  {
    my ($self,$x)=@_;
    $self->{begin}=$x;
  }
#---------------------------------------------------------------------
#   $transcript->setEnd($x);
sub setEnd
  {
    my ($self,$x)=@_;
    $self->{end}=$x;
  }
#---------------------------------------------------------------------
#   $genomicCoord=$transcript->mapToGenome($transcriptCoord);
sub mapToGenome
  {
    my ($self,$transcriptCoord)=@_;
    my $original=$transcriptCoord;
    my $numExons=$self->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $exon=$self->getIthExon($i);
      my $exonLen=$exon->getLength();
      if($transcriptCoord<$exonLen) {
	return $self->getStrand() eq "+" ?
	  $exon->getBegin()+$transcriptCoord :
	    $exon->getEnd()-$transcriptCoord-1;
      }
      $transcriptCoord-=$exonLen;
    }
    my $id=$self->getID();
    die "coordinate is beyond transcript end: $original ($id)";
  }
#---------------------------------------------------------------------
#   $transcriptCoord=$transcript->mapToTranscript($genomicCoord);
sub mapToTranscript
  {
    my ($self,$genomicCoord)=@_;
    my $numExons=$self->numExons();
    my $transcriptCoord=0;
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $exon=$self->getIthExon($i);
      if($exon->containsCoordinate($genomicCoord)) {
	return $self->getStrand() eq "+" ?
	  $transcriptCoord+$genomicCoord-$exon->getBegin() :
	    $transcriptCoord+$exon->getEnd()-$genomicCoord-1;
	}
      my $exonLen=$exon->getLength();
      $transcriptCoord+=$exonLen;
    }
    my $id=$self->getID();
    die "transcript $id does not contain genomic coordinate $genomicCoord";
  }
#---------------------------------------------------------------------
#   $exon=$transcript->getExonContaining($genomicCoord);
sub getExonContaining
  {
    my ($self,$genomicCoord)=@_;
    my $numExons=$self->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $exon=$self->getIthExon($i);
      if($exon->containsCoordinate($genomicCoord)) { return $exon }
    }
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#   $copy=$transcript->copy();
sub copy 
  {
    my ($self)=@_;
    my $new=new Transcript;
    %$new=%$self;
    $new->{exons}=[];
    my $numExons=$self->numExons();
    for(my $i=0 ; $i<$numExons ; ++$i) {
      my $newExon=$self->getIthExon($i)->copy();
      $newExon->{transcript}=$new;
      $new->addExon($newExon);
    }
    return $new;
  }
#---------------------------------------------------------------------
#   $transcript->reverseExonOrder();
sub reverseExonOrder {
  my ($self)=@_;
  my $exons=$self->{exons};
  my $a=[];
  my $numExons=@$exons;
  for(my $i=$numExons-1 ; $i>=0 ; --$i) {
    push @$a,$exons->[$i];
  }
  $self->{exons}=$a;
}
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

1;

