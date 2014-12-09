package PromerMatch;
use strict;

######################################################################
#
# PromerMatch.pm bmajorostigr.org 8/22/2003
#
# 
# 
#
# Attributes:
#   begin1
#   end1
#   begin2
#   end2
#   length1
#   length2
#   strand1
#   strand2
#   percentIdentity
#   percentSimilarity
#   percentStops
#   substrate1
#   substrate2
#   frame1
#   frame2
# Methods:
#   $promerMatch=new PromerMatch($begin1,$end1,$begin2,$end2,$length1,$length2,
#                                $strand1,$strand2,$identity,$similarity,
#                                $stops,$subst1,$subst2,$frame1,$frame2);
#   $orientation=$promerMatch->getOrientation(); # '+' or '-'
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $promerMatch=new PromerMatch($begin1,$end1,$begin2,$end2,$length1,$length2,
#                                $strand1,$strand2,$identity,$similarity,
#                                $stops,$subst1,$subst2,$frame1,$frame2);
sub new
{
  my ($class,$begin1,$end1,$begin2,$end2,$length1,$length2,$strand1,
      $strand2,$identity,$sim,$stops,$subst1,$subst2,$frame1,$frame2)=@_;
  
  my $self=
    {
     begin1=>$begin1,
     end1=>$end1,
     begin2=>$begin2,
     end2=>$end2,
     length1=>$length1,
     length2=>$length2,
     strand1=>$strand1,
     strand2=>$strand2,
     percentIdentity=>$identity,
     percentSimilarity=>$sim,
     percentStops=>$stops,
     substrate1=>$subst1,
     substrate2=>$subst2,
     frame1=>$frame1,
     frame2=>$frame2
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $orientation=$promerMatch->getOrientation(); # '+' or '-'
sub getOrientation
  {
    my ($self)=@_;
    if(($self->{end1}-$self->{begin1})/
       ($self->{end2}-$self->{begin2}) > 0) {return '+'}
    return '-';
  }
#---------------------------------------------------------------------
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

