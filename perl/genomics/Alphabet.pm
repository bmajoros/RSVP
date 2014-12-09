package Alphabet;
use strict;

######################################################################
#
# Alphabet.pm bmajorostigr.org 10/10/2003
#
# 
# 
#
# Attributes:
#   @members
# Methods:
#   $alphabet=new Alphabet("ATCG");
#   my $int=$alphabet->symbolToIndex($symbol);
#   my $symbol=$alphabet->indexToSymbol($int);
#   my $size=$alphabet->size();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $alphabet=new Alphabet("ATCG");
sub new
{
  my ($class,$initString)=@_;

  my @members=split//,$initString;
  my $self={members=>\@members};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   my $int=$alphabet->symbolToIndex($symbol);
sub symbolToIndex
  {
    my ($self,$symbol)=@_;
    my $members=$self->{members};
    my $size=@$members;
    for(my $i=0 ; $i<$size ; ++$i)
      {
	my $member=$members->[$i];
	if($member eq $symbol) {return $i}
      }
    die "symbol $symbol not found in symbolToIndex()";
  }
#---------------------------------------------------------------------
#   my $symbol=$alphabet->indexToSymbol($int);
sub indexToSymbol
  {
    my ($self,$index)=@_;
    return $self->{members}->[$index];
  }
#---------------------------------------------------------------------
#   my $size=$alphabet->size();
sub size
  {
    my ($self)=@_;
    return 0+@{$self->{members}};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

