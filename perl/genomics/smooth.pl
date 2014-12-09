#!/usr/bin/perl
use strict;

die "smooth.pl <window-size> <num-iterations>\n" unless @ARGV==2;
my ($windowSize,$numIterations)=@ARGV;

my @histogram;
while(<STDIN>)
  {
    if(/(\S+)\s+(\S+)/)
      {
	next if($2>1000 || $2<-1000);
	push @histogram,[$1,$2];
      }
  }
for(my $i=0 ; $i<$numIterations ; ++$i)
  {
    @histogram=@{smooth($windowSize,\@histogram)};
  }
my $L=@histogram;
for(my $i=0 ; $i<$L ; ++$i)
  {
    my $pair=$histogram[$i];
    my ($x,$y)=@$pair;
    print "$x\t$y\n";
  }

my $FALSE=0;
my $TRUE=1;
my $NUM_LEFT_SKIP_BUCKETS=1;
my $SMOOTH_LEFT_MARGIN=$FALSE;

sub smooth
  {
    my ($windowSize,$histogram)=@_;
    my $n=@$histogram;
    my $halfWindow=int($windowSize/2);
    my $otherHalf=$windowSize-$halfWindow; # in case it's an odd number
    my $first=$halfWindow;
    my $last=$n-1-$otherHalf;
    my $newHistogram=[];

    # Handle the leftmost bins (too close to edge to use a full window)
    my $boundarySum;
    if($SMOOTH_LEFT_MARGIN)
      {
	for(my $i=0 ; $i<$windowSize ; ++$i) 
	  {
	    my $pair=$histogram->[$i];
	    my $y=(defined($pair) ? $pair->[1] : 0);
	    $boundarySum+=$y;
	  }
	my $boundaryAve=$boundarySum/$windowSize;
	for(my $i=0 ; $i<$first ; ++$i) 
	  {
	    $newHistogram->[$i]=$histogram->[$i];
	    $newHistogram->[$i]->[1]=$boundaryAve;
	  }
      }
    else
      {
	for(my $i=0 ; $i<$first ; ++$i) 
	  {
	    $newHistogram->[$i]=$histogram->[$i];
	  }
      }

    # Handle the rightmost bins (too close to edge to use a full window)
    $boundarySum=0;
    for(my $i=$last+1 ; $i<$n ; ++$i) 
      {
	my $pair=$histogram->[$i];
	my $y=(defined($pair) ? $pair->[1] : 0);
	$boundarySum+=$y;
      }
    my $boundaryAve=$boundarySum/$windowSize;
    for(my $i=$last+1 ; $i<$n ; ++$i) 
      {
	$newHistogram->[$i]=$histogram->[$i];
	$newHistogram->[$i]->[1]=$boundaryAve;
      }

    # Handle the main part of the histogram
    for(my $i=$first ; $i<=$last ; ++$i)
      {
	my $pair=$histogram->[$i];
	my ($x,$y)=@$pair;
	for(my $j=0 ; $j<$halfWindow ; ++$j)
	  {
	    my $pair=$histogram->[$i-1-$j];
	    my ($leftX,$leftY)=(defined($pair) ? @$pair : (0,0));
	    $y+=$leftY;
	  }
	for(my $j=0 ; $j<$otherHalf ; ++$j)
	  {
	    my $pair=$histogram->[$i+1+$j];
	    my ($rightX,$rightY)=(defined($pair) ? @$pair : (0,0));
	    $y+=$rightY;
	  }
	$y/=($windowSize+1);
	$newHistogram->[$i]=[$x,$y];
      }
    return $newHistogram;
  }
