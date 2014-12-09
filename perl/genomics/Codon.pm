package Codon;
use strict;
use Carp;

######################################################################
#
# Codon.pm 
# bmajoros@tigr.org 7/8/2002
#
# Represents a single trinucleotide at a certain place within an exon.
# Stores its location both relative to the transcript and to the
# genomic axis.  All coordinates are zero-based and space-based.  On
# the forward strand, the absoluteCoord represents the position of the
# leftmost base in the codon, whereas on the reverse strand,
# the absoluteCoord points to the base immediately following the third
# base of the codon.  Thus, on the minus strand, a substring operation
# on the genomic axis should begin at absoluteCoord-3.  RelativeCoords
# are much simpler: a substring operation on the transcript always
# uses the relativeCoord as-is, regardless of strand.
#
# Attributes:
#   string triplet
#   int absoluteCoord : relative to genomic axis
#   int relativeCoord : relative to current exon
#   bool isInterrupted : exon ends before codon is complete
#   int basesInExon :  how many bases of this codon are in this exon (1-3)
#   Exon *exon : which exon contains this codon
# Methods:
#   $codon=new Codon($exon,$triplet,$relative,$absolute,$isInterrupted);
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $codon=new Codon($exon,$triplet,$relative,$absolute,$isInterrupted);
sub new
{
  my ($class,$exon,$triplet,$relative,$absolute,$isInterrupted)=@_;
  confess unless defined $relative;

  my $self=
    {
     triplet=>$triplet,
     exon=>$exon,
     relativeCoord=>$relative,
     absoluteCoord=>$absolute,
     isInterrupted=>$isInterrupted,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

