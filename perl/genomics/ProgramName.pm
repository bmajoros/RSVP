package ProgramName;
use strict;

######################################################################
#
# ProgramName.pm bmajorostigr.org 3/5/2005
#
# Uses $0 to get the program name and then strips the path from the
# filename so that the usage statement won't be long on account of
# the program existing in some deep directory.
#
# Methods:
#   $name=ProgramName::get();
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
#   $name=ProgramName::get();
sub get
  {
    $0=~/([^\/]+)\s*$/;
    return $1;
  }





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

