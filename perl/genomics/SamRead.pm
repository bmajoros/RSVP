package SamRead;
use strict;

######################################################################
#
# SamRead.pm bmajoros@duke.edu 7/11/2011
#
# 
# 
#
# Attributes:
#   id : read identifier
#   flags
#   substrate : identifier of scaffold/chromosome
#   pos : 0-based coordinate of read on substrate
#   strand : '+' or '-'
#   mqual
#   cigar : CIGAR string (alignment)
#   rnext
#   pnext
#   tlen
#   seq : read sequence
#   qual : sequence quality string
# Methods:
#   $samRead=new SamRead($id,$flags,$substrate,$pos,$mqual,$cigar,
#                        $rnext,$pnext,$tlen,$seq,$qual);
#   $strand=$read->getStrand();
#   $id=$read->getID();
#   $substrate=$read->getSubstrate();
#   $pos=$read->getBegin();  # 0-based inclusive
#   $endPos=$read->getEnd(); # 0-based non-inclusive (after the last base)
#   $cigar=$read->getCigarString();
#   $seq=$read->getSeq();
#   $qual=$read->getQualityString();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class,$id,$flags,$substrate,$pos,$mqual,$cigar,$rnext,$pnext,$tlen,$seq,$qual)=@_;
  my $strand=($flags & 0x10) ? '-' : '+';
  my $self=
    {
     class=>$class,
     id=>$id,
     strand=>$strand,
     flags=>$flags,
     substrate=>$substrate,
     pos=>$pos,
     mqual=>$mqual,
     cigar=>$cigar,
     rnext=>$rnext,
     pnext=>$pnext,
     tlen=>$tlen,
     seq=>$seq,
     qual=>$qual
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $strand=$read->getStrand();
sub getStrand
  {
    my ($self)=@_;
    return $self->{strand};
  }
#---------------------------------------------------------------------
#   $id=$read->getID();
sub getID()
  {
    my ($self)=@_;
    return $self->{id};
  }
#---------------------------------------------------------------------
#   $substrate=$read->getSubstrate();
sub getSubstrate
  {
    my ($self)=@_;
    return $self->{substrate};
  }
#---------------------------------------------------------------------
#   $pos=$read->getBegin();
sub getBegin
  {
    my ($self)=@_;
    return $self->{pos};
  }
#---------------------------------------------------------------------
#   $endPos=$read->getEnd();
sub getEnd
  {
    my ($self)=@_;
    my $begin=$self->getBegin();
    my $cigar=$self->getCigarString();
    #42M105N8M
    my @lengths=split/[MN]/,$cigar;
    my $L=0;
    foreach my $len (@lengths) {$L+=$len}
    return $begin+$L;

  }
#---------------------------------------------------------------------
#   $cigar=$read->getCigarString();
sub getCigarString
  {
    my ($self)=@_;
    return $self->{cigar};
  }
#---------------------------------------------------------------------
#   $seq=$read->getSeq();
sub getSeq
  {
    my ($self)=@_;
    return $self->{seq};
  }
#---------------------------------------------------------------------
#   $qual=$read->getQualityString();
sub getQualityString
  {
    my ($self)=@_;
    return $self->{qual};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

