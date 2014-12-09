package LinearRegressor;
use strict;

######################################################################
#
# LinearRegressor.pm bmajoros
#
# 
# 
#
# Attributes:
#
# Methods:
#   $regressor=new LinearRegressor();
#   ($slope,$intercept,$r,$r2)=$regressor->regress(\@x,\@y);
#      # $r=correlation coeff., $r2=coeff. of determination
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
#   ($slope,$intercept,$r,$r2)=$regressor->regress(\@x,\@y);
sub regress
  {
    my ($self,$xs,$ys)=@_;

    my ($sumX,$sumXX,$sumY,$sumXY,$sumYY);

    my $n=@$xs;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $x=$xs->[$i];
	my $y=$ys->[$i];
	
	$sumX+=$x;
	$sumXX+=$x*$x;
	$sumY+=$y;
	$sumXY+=$x*$y;
	$sumYY+=$y*$y;
      }
    
    my $Sxx=$sumXX-$sumX*$sumX/$n;
    my $Sxy=$sumXY-$sumX*$sumY/$n;
    my $Syy=$sumYY-$sumY*$sumY/$n;
    
    my $xBar=$sumX/$n;
    my $yBar=$sumY/$n;
    
    if($Sxx==0.0) { die "Sxx=zero" }
    my $slope=$Sxy/$Sxx;
    my $intercept=$yBar-$slope*$xBar;
    
    my $r=($Sxx*$Syy!=0 ? $Sxy/sqrt($Sxx*$Syy) : 1);
    my $r2=$r*$r;
    
    return ($slope,$intercept,$r,$r2);
  }





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

