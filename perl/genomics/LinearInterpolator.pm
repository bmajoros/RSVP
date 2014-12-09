package LinearInterpolator;
use strict;

######################################################################
#
# LinearInterpolator.pm bmajoros
#
# 
# 
#
# Attributes:
#
# Methods:
#   $interp=new LinearInterpolator();
#   my $y=$interp->interpolate($x1,$y1,$x2,$y2,$x);
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self={};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   my $y=$interp->interpolate($x1,$y1,$x2,$y2,$x);
sub interpolate
  {
    my ($self,$x1,$y1,$x2,$y2,$x)=@_;
    my $dy=$y2-$y1;
    my $dx=$x2-$x1;
    my $m=$dy/$dx;
    my $b=$y1;
    my $y=$m*($x-$x1)+$b;
    return $y;
  }





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

