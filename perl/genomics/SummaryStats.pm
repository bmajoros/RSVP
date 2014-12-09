package SummaryStats;
use strict;
use Carp;

######################################################################
#
# SummaryStats.pm bmajoros 2/19/2001
#
# 
# 
# ($mean,$stddev,$min,$max)=SummaryStats::summaryStats(\@array);
# ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@array);
# $sum=SummaryStats::sum(\@array);
# $r=SummaryStats::correlation(\@array1,\@array2);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
# ($mean,$stddev,$min,$max)=SummaryStats::summaryStats(\@array);
sub summaryStats
  {
    my ($array)=@_;
    my $n=@$array;
    my ($minX,$maxX);
    my $sumX=0;
    my $sumXX=0;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $x=$array->[$i];
	$sumX+=$x;
	$sumXX+=($x*$x);
	$minX=$maxX=$x if $i==0;
	$minX=$x if $x<$minX;
	$maxX=$x if $x>$maxX;
      }
    my $meanX=$sumX/$n;
    #confess "n\=$n in summaryStats()\n" unless $n>1;
    my $varX=$n>1 ? ($sumXX-$sumX*$sumX/$n)/($n-1) : undef;
    if($varX<0) {$varX=0}
    my $stddevX=sqrt($varX);
    return ($meanX,$stddevX,$minX,$maxX);
  }
#---------------------------------------------------------------------
# $sum=SummaryStats::sum(\@array);
sub sum  {
  my ($array)=@_;
  my $s=0.0;
  my $n=@$array;
  for(my $i=0 ; $i<$n ; ++$i) {
    $s+=$array->[$i];
  }
  return $s;
}
#---------------------------------------------------------------------
# $r=SummaryStats::correlation(\@array1,\@array2);
sub correlation
  {
    my ($Xs,$Ys)=@_;

    my $sumX=0.0;
    my $sumY=0.0;
    my $sumXY=0.0;
    my $sumXX=0.0;
    my $sumYY=0.0;
    my $n=@$Xs;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $x=$Xs->[$i];
	my $y=$Ys->[$i];

	$sumX+=$x;
	$sumY+=$y;
	$sumXY+=($x*$y);
	$sumXX+=($x*$x);
	$sumYY+=($y*$y);
      }
    
    my $r=($sumXY-$sumX*$sumY/$n)/
      sqrt(($sumXX-$sumX*$sumX/$n)*($sumYY-$sumY*$sumY/$n));
    return $r;
  }
#---------------------------------------------------------------------
# ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@array);
sub roundedSummaryStats
  {
    my ($array)=@_;
    my ($mean,$stddev,$min,$max)=SummaryStats::summaryStats($array);
    $mean=int(100*$mean+5/9)/100;
    $stddev=int(100*$stddev+5/9)/100;
    $min=int(100*$min+5/9)/100;
    $max=int(100*$max+5/9)/100;
    return ($mean,$stddev,$min,$max);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

