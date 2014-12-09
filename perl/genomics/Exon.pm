package Exon;
use strict;
use Feature;

######################################################################
#
# Exon.pm
#
# bmajoros@tigr.org 7/8/2002
#
# A pair of coordinates representing the location of an exon on a
# genomic axis.  The first coordinate (begin) is always less than
# the second coordinate (end), regardless of strand.  The sequence
# is only loaded if loadExonSequences() or loadTranscriptSeq() is sent
# to the parent transcript.  The sequence is automatically reverse-
# complemented if on the reverse strand.  Has an order index telling
# you which exon in the transcript this is, from the first (0) to
# the last, in translation order.  Thus, the first exon (order=0
# is actually rightmost when on the reverse strand.  Note that frame
# (actually phase) denotes the phase of the first base of an exon, but
# on the reverse strand, the first base is actually the rightmost.
#
# Attributes:
#   begin : the left end of the exon (zero-based)
#   end   : one base past the right end of the exon (zero-based)
#   order : which exon this is (0,1,2,...), in translation order
#   sequence : NOT LOADED BY DEFAULT! (but automatically reverse-
#              complemented when it is loaded)
#   transcript : the parent transcript
#   frame : frame of first base in exon
#   type : type of object, usually "exon" or "internal-exon", etc...
#   score
#   strand
#   substrate
# Methods:
#   $exon=new Exon($begin,$end,$transcript);
#   $exon->containsCoordinate($x) : boolean
#   $length=$exon->getLength();
#   $exon->reverseComplement($seqLen);
#   $exon->trimInitialPortion($numBases);
#   $exon->trimFinalPortion($numBases);
#   $strand=$exon->getStrand(); # "+" or "-"
#   $bool=$exon->overlaps($otherExon)
#   $sequence=$exon->getSequence();
#   $transcript=$exon->getTranscript();
#   $gff=$exon->toGff();
#   $begin=$exon->getBegin();
#   $end=$exon->getEnd();
#   $frame=$exon->getFrame();
#   $exon->setFrame($frame);
#   $type=$exon->getType();
#   $exon->setType($type);
#   $exon->setScore($score);
#   $score=$exon->getScore();
#   $exon->shiftCoords($delta);
#   $exon->setStrand($strand);
#   $substrate=$exon->getSubstrate();
#   $exon->setSubstrate($substrate);
#   $exon->setBegin($begin);
#   $exon->setEnd($end);
######################################################################

#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $exon=new Exon($begin,$end,$transcript);
sub new
{
  my ($class,$begin,$end,$transcript)=@_;
  
  my $self=
    {
     begin=>$begin,
     end=>$end,
     transcript=>$transcript,
     score=>".",
     strand=>$transcript->getStrand(),
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $exon->setBegin($begin);
sub setBegin
  {
    my ($self,$begin)=@_;
    $self->{begin}=$begin;
  }
#---------------------------------------------------------------------
#   $exon->setEnd($end);
sub setEnd
  {
    my ($self,$end)=@_;
    $self->{end}=$end;
  }
#---------------------------------------------------------------------
#   $exon->containsCoordinate($x) : boolean
sub containsCoordinate
  {
    my ($self,$x)=@_;
    return $x>=$self->{begin} && $x<$self->{end};
  }
#---------------------------------------------------------------------
#   $length=$exon->getLength();
sub getLength
  {
    my ($self)=@_;
    return $self->{end}-$self->{begin};
  }
#---------------------------------------------------------------------
#   $exon->trimInitialPortion($numBases);
#
# Trims a certain number of bases from the translationally early
# part of the exon.
sub trimInitialPortion
  {
    my ($self,$numBases)=@_;
    if($self->getStrand() eq "+")
      {
	$self->{begin}+=$numBases;
      }
    else # strand eq "-"
      {
	$self->{end}-=$numBases;
      }
    my $sequence=$self->{sequence};
    if(defined $sequence)
      {
	$self->{sequence}=
	  substr($sequence,$numBases,length($sequence)-$numBases);
      }
  }
#---------------------------------------------------------------------
#   $exon->trimFinalPortion($numBases);
#
# Trims a certain number of bases from the translationally late
# part of the exon.
sub trimFinalPortion
  {
    my ($self,$numBases)=@_;
    if($self->getStrand() eq "+")
      {
	$self->{end}-=$numBases;
      }
    else # strand eq "-"
      {
	$self->{begin}+=$numBases;
      }
    my $sequence=$self->{sequence};
    if(defined $sequence)
      {
	$self->{sequence}=
	  substr($sequence,0,length($sequence)-$numBases);
      }
  }
#---------------------------------------------------------------------
#   $strand=$exon->getStrand(); # "+" or "-"
sub getStrand
  {
    my ($self)=@_;
    return $self->{transcript}->{strand};
  }
#---------------------------------------------------------------------
#   $bool=$exon->overlaps($otherExon)
sub overlaps
  {
    my($this,$otherExon)=@_;
    return 
      $this->{begin}<$otherExon->{end} &&
      $otherExon->{begin}<$this->{end};
  }
#---------------------------------------------------------------------
#   $sequence=$exon->getSequence(); 
sub getSequence
  {
    my ($self)=@_;
    return $self->{sequence};
  }
#---------------------------------------------------------------------
#   $transcript=$exon->getTranscript();
sub getTranscript
  {
    my ($self)=@_;
    return $self->{transcript};
  }
#---------------------------------------------------------------------
#   $begin=$exon->getBegin();
sub getBegin
  {
    my ($self)=@_;
    return $self->{begin};
  }
#---------------------------------------------------------------------
#   $end=$exon->getEnd();
sub getEnd
  {
    my ($self)=@_;
    return $self->{end};
  }
#---------------------------------------------------------------------
#   $frame=$exon->getFrame();
sub getFrame
  {
    my ($self)=@_;
    my $frame=$self->{frame};
    return (defined($frame) ? $frame : ".");
  }
#---------------------------------------------------------------------
#   $type=$exon->getType();
sub getType
  {
    my ($self)=@_;
    return $self->{type};
  }
#---------------------------------------------------------------------
#   $gff=$exon->toGff();
sub toGff
  {
    my ($self)=@_;
    my $transcript=$self->getTranscript();
    my $substrate=$self->getSubstrate();
    my $begin=$self->getBegin()+1; # convert to 1-based coordinate system (1/B)
    my $end=$self->getEnd();
    my $type=$self->getType() || ".";
    my $source=($transcript ? $transcript->getSource() || "." : ".");
    my $strand=$self->getStrand();
    my $frame=$self->getFrame();
    my $score=$self->getScore();
    my $transcriptId=($transcript ? $transcript->getID() : ".");
    my $geneId=($transcript ? $transcript->getGeneId() : ".");
    my $gff="$substrate\t$source\t$type\t$begin\t$end\t$score\t$strand\t$frame\ttranscript_id \"$transcriptId\"; gene_id \"$geneId\";\n";
    return $gff;
  }
#---------------------------------------------------------------------
#   $exon->setScore($score);
sub setScore
  {
    my ($self,$score)=@_;
    $self->{score}=$score;
  }
#---------------------------------------------------------------------
#   my $score=$exon->getScore();
sub getScore
  {
    my ($self)=@_;
    my $score=$self->{score};
    return (defined($score) ? $score : ".");
  }
#---------------------------------------------------------------------
#   $exon->setType($type);
sub setType
  {
    my ($self,$type)=@_;
    $self->{type}=$type;
  }
#---------------------------------------------------------------------
#   $exon->shiftCoords($delta);
sub shiftCoords
  {
    my ($self,$delta)=@_;
    $self->{begin}+=$delta;
    $self->{end}+=$delta;
  }
#---------------------------------------------------------------------
#   $exon->setFrame($frame);
sub setFrame
  {
    my ($self,$frame)=@_;
    $self->{frame}=$frame;
  }
#---------------------------------------------------------------------
#   $strand=$exon->getStrand();
sub getStrand
  {
    my ($self)=@_;
    my $strand=$self->{strand};
    if(defined $strand) {return $strand}
    return $self->{transcript}->getStrand();
  }
#---------------------------------------------------------------------
#   $exon->setStrand($strand);
sub setStrand
  {
    my ($self,$strand)=@_;
    $self->{strand}=$strand;
  }
#---------------------------------------------------------------------
#   $substrate=$exon->getSubstrate();
sub getSubstrate
  {
    my ($self)=@_;
    my $substrate=$self->{substrate};
    if(defined $substrate) {return $substrate}
    return $self->{transcript}->getSubstrate();
  }
#---------------------------------------------------------------------
#   $exon->setSubstrate($substrate);
sub setSubstrate
  {
    my ($self,$substrate)=@_;
    $self->{substrate}=$substrate;
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
#   $exon->reverseComplement($seqLen);
sub reverseComplement
  {
    my ($self,$seqLen)=@_;
    my $begin=$self->getBegin();
    my $end=$self->getEnd();
    $self->{begin}=$seqLen-$end;
    $self->{end}=$seqLen-$begin;
    $self->{strand}=compStrand($self->{strand});
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

