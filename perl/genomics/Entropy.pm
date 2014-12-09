package Entropy;
use strict;

######################################################################
#
# Entropy.pm bmajoros 7/19/2002
#
# 
# 
#
# Attributes:
#
# Methods:
#   $entropy=Entropy::entropy($sequence); # pass pointer to sequence!
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub entropy
  {
    my ($seq)=@_;
    my $len=length($$seq);
    my (%counts,$n);
    for(my $i=0 ; $i<$len ; ++$i)
      {
	my $residue=substr($$seq,$i,1);
	++$n;
	++$counts{$residue};
      }
    my @alphabet=keys %counts;
    my $H=0;
    foreach my $symbol (@alphabet)
      {
	my $p=$counts{$symbol}/$n;
	next unless $p>0;
	$H-=$p*ln($p);
      }
    return $H;
  }
#---------------------------------------------------------------------
sub ln
  {
    my ($x)=@_;
    return log($x)/log(2);
  }





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

