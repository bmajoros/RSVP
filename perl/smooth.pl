#!/usr/bin/perl
#Bill's code

use strict;

my $ITERATIONS=1000;
my $WINDOW_SIZE=30;
my $SMOOTH_LEFT_MARGIN=1;

my @hist;
while(<STDIN>) {
 chomp;
 next unless(/(\d+)\s+(\S+)$/);
 push @hist,[$1,$2];
}
my $smoothed=\@hist;
for(my $i=0 ; $i<$ITERATIONS ; ++$i) {
 $smoothed=smooth($WINDOW_SIZE,$smoothed);
}
my $L=@$smoothed;
for(my $i=0 ; $i<$L ; ++$i) {
 my ($x,$y)=@{$smoothed->[$i]};
 print "$x\t$y\n";
}




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
         {
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
         }
       $newHistogram->[$i]=[$x,$y];
     }

   return $newHistogram;
 }
