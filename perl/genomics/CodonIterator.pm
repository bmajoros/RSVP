package CodonIterator;
use strict;
use Codon;
use Exon;
use Transcript;
use Translation;

######################################################################
#
# CodonIterator.pm 
# bmajoros@tigr.org 7/8/2002
#
# Iterates through the codons of a transcript.  Each codon that it
# produces will have an absolute and relative coordinate, an
# indicator of which exon contains it, and whether the codon is
# interrupted by the end of the exon.
#
# Attributes:
#   transcript : Transcript
#   exon : which exon we are processing
#   relative : current position relative to current exon's 5-prime end
#   absolute : current position relative to genomic axis
#   stopCodons : pointer to hash containing stop codon sequences
# Methods:
#   $iterator=new CodonIterator($transcript,\$axisSequence,\%stopCodons);
#   $codon=$iterator->nextCodon(); # or undef if no more
#   $codons=$iterator->getAllCodons();
#
# Private methods:
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $iterator=new CodonIterator($transcript,\$axisSequence,\%stopCodons);
sub new
{
  my ($class,$transcript,$axisSequenceRef,$stopCodons)=@_;

  my $self=
    {
     transcript=>$transcript,
     stopCodons=>$stopCodons
    };
  bless $self,$class;

  # Advance to the exon containing the start codon
  my $exons=$transcript->{exons};
  my $strand=$transcript->{strand};
  my $startCodon=$transcript->{startCodonAbsolute};
  if(defined($startCodon))
    {
      die unless length($$axisSequenceRef)>0; ###
      $transcript->loadExonSequences($axisSequenceRef);
      my $numExons=@$exons;
      my $exon;
      for(my $i=0 ; $i<$numExons ; ++$i)
	{
	  $exon=$exons->[$i];
	  if($exon->containsCoordinate($startCodon) ||
	     $exon->containsCoordinate($startCodon-1)) {last}
	}
      die "$startCodon not found" unless $exon->containsCoordinate($startCodon)
	|| $exon->containsCoordinate($startCodon-1);
      $self->{exon}=$exon;
      $self->{relative}=
	($strand eq "+" ?
	 $startCodon-$exon->{begin} :
	 $exon->{end}-$startCodon);
      $self->{absolute}=$startCodon;
    }
  else
    {
      my $exon=$exons->[0];
      $self->{exon}=$exon;
      my $frame=$exon->getFrame();
      die unless defined $frame;
      my $add=(3-$frame)%3;
      if($strand eq "+")
	{
	  $self->{absolute}=$exon->getBegin()+$add;
	  $self->{relative}=$add;
	}
      else
	{
	  $self->{absolute}=$exon->getEnd()-1-$add;
	  $self->{relative}=$add;
	}
    }

  return $self;
}
#---------------------------------------------------------------------
#   $codons=$iterator->getAllCodons();
sub getAllCodons
  {
    my ($self)=@_;
    my @codons;
    while(1)
      {
	my $codon=$self->nextCodon();
	last unless defined $codon;
	push @codons,$codon;
      }
    return \@codons;
  }
#---------------------------------------------------------------------
#   $codon=$iterator->nextCodon();
sub nextCodon
  {
    my ($self)=@_;
    my $isStopCodon=$self->{stopCodons};
    my $exon=$self->{exon};
    if(!defined($exon)) {return undef}
    my $exonSeq=$exon->{sequence};
    my $relative=$self->{relative};
    my $absolute=$self->{absolute};
    my $exonBegin=$exon->{begin};
    my $exonEnd=$exon->{end};
    my $exonLen=$exonEnd-$exonBegin;
    my $transcript=$self->{transcript};
    my $transcriptId=$transcript->getID();
    my $strand=$transcript->{strand};
    my $isInterrupted=$relative>$exonLen-3;
    my $triplet;
    my $thisExonContrib=3;
    if($isInterrupted)
      {
	$thisExonContrib=$exonLen-$relative;
	my $nextExonContrib=3-$thisExonContrib;
	$triplet=substr($exonSeq,$relative,$thisExonContrib);
	my $transcript=$self->{transcript};
	my $thisExonOrder=$exon->{order};
	die "exon has no order" unless defined $thisExonOrder;
	my $exons=$transcript->{exons};
	my $numExons=@$exons;
	my $nextExonOrder=$thisExonOrder+1;
	if($nextExonOrder>=$numExons)
	  {undef $self->{exon};return undef} # end of iteration
	else
	  {
	    my $nextExon=$exons->[$nextExonOrder];
	    $exonSeq=$nextExon->{sequence};
	    $triplet.=substr($exonSeq,0,$nextExonContrib);
	    if(length($triplet)!=3)### debugging
	      {
		my $realSeqLen1=length $exon->{sequence};
		my $exonLen1=$exon->getLength();
		my $realSeqLen2=length $nextExon->{sequence};
		my $exonLen2=$nextExon->getLength();

		die "Error in transcript $transcriptId: nextContrib=$nextExonContrib triplet=\"$triplet\" exonLen1=$exonLen1 realLen1=$realSeqLen1 exonLen2=$exonLen2 realLen2=$realSeqLen2";
	      }
	    $self->{relative}=$nextExonContrib;
	    $self->{absolute}=
	      ($strand eq "+" ? 
	       $nextExon->{begin}+$nextExonContrib :
	       $nextExon->{end}-$nextExonContrib); ###-1);
	    $self->{exon}=$nextExon;
	  }
      }
    else # codon was not interrupted by end of exon
      {
	$triplet=substr($exonSeq,$relative,3);
	if($isStopCodon->{$triplet})
	  {undef $self->{exon}} # end of iteration
	else
	  {
	    $self->{relative}+=3;
	    $self->{absolute}+=($strand eq "+" ? 3 : -3);
	    if($self->{relative}>=$exon->{end}-$exon->{begin})
	      {
		my $transcript=$self->{transcript};
		my $thisExonOrder=$exon->{order};
		my $exons=$transcript->{exons};
		my $numExons=@$exons;
		my $nextExonOrder=$thisExonOrder+1;
		if($nextExonOrder>=$numExons)
		  {undef $self->{exon}} # end of iteration
		else
		  {
		    my $nextExon=$exons->[$nextExonOrder];
		    $self->{exon}=$nextExon;
		    $self->{relative}=0;
		    $self->{absolute}=
		      ($strand eq "+" ? 
		       $nextExon->{begin} : 
		       $nextExon->{end});
		  }
	      }
	  }
      }
    my $codon=new Codon($exon,$triplet,$relative,$absolute,$isInterrupted);
    $codon->{basesInExon}=$thisExonContrib;
    return $codon;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

1;

