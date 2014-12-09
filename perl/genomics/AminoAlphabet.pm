package AminoAlphabet;
use strict;
use Alphabet;
use vars qw(@ISA);

@ISA=qw(Alphabet);

######################################################################
#
# AminoAlphabet.pm bmajorostigr.org 10/10/2003
#
# 
# 
#
# Attributes:
#
# Methods:
#   $aminoAlphabet=new AminoAlphabet();
#
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;

  my $self=new Alphabet("ARNDCQEGHILKMFPSTWYVBZX*");
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

