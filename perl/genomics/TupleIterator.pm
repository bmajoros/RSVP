package TupleIterator;
use strict;

######################################################################
#
# TupleIterator.pm bmajoros
#
# Enumerates all N-tuples for a given N where each element is less
# than all those to its right.  For example, (3,4,7) is a valid
# 3-tuple, but (1,1,5) and (4,3,8) are not.  $N in the constructor
# is the length of the tuple, and $maxValue constrains the values
# in each position (which are nonnegative integers).
#
# Attributes:
#   N : int
#   maxValue : int
#   tuple : array reference
# Methods:
#   $iterator=new TupleIterator($N,$maxValue);
#   $arrayRef=$iterator->currentTuple(); # returns undef when finished
#   $iterator->advance(); # returns undef when finished
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $iterator=new TupleIterator($N,$maxValue);
sub new
{
  my ($class,$N,$maxValue)=@_;

  my $tuple=[];
  for(my $i=0 ; $i<$N ; ++$i)
    {
      $tuple->[$i]=$i;
    }

  my $self=
    {
     N => $N,
     maxValue => $maxValue,
     tuple => $tuple,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $arrayRef=$iterator->currentTuple();
sub currentTuple
  {
    my ($self)=@_;

    return $self->{tuple};
  }
#---------------------------------------------------------------------
#   $iterator->advance();
sub advance
  {
    my ($self)=@_;

    my $N=$self->{N};
    $self->inc($N-1);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
sub inc
  {
    my ($self,$m)=@_;

    if($m<0)
      {
	undef $self->{tuple};
	return undef;
      }

    my $maxInt=$self->{maxInt};
    my $tuple=$self->{tuple};
    my $maxValue=$self->{maxValue};
    my $N=$self->{N};
    
    ++$tuple->[$m];
    my $Mmax=$maxValue+$m-$N+1;
    if($tuple->[$m]>$Mmax)
      {
	return $self->inc($m-1);
      }
    else
      {
	my $x=$tuple->[$m]+1;
	for(my $i=$m+1 ; $i<$N ; ++$i)
	  {
	    $tuple->[$i]=$x;
	    ++$x;
	  }
      }
    return 1;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


1;

