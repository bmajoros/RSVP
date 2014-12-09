package NgramIterator;
use strict;

######################################################################
#
# NgramIterator.pm bmajoros
#
# 
# 
#
# Attributes:
#   array ngram : array of integer indices into alphabet string
#   string alphabet
# Methods:
#   $ngramIterator=new NgramIterator("ATCG",$N);
#   $string=$ngramIterator->nextString(); # returns undef if no more
#   
# Private methods:
#   $self->ngramToString();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $ngramIterator=new NgramIterator("ATCG",$N);
sub new
{
  my ($class,$alphabet,$N)=@_;

  my $ngram=[];
  my $alphaSize=length $alphabet;
  for(my $i=0 ; $i<$N-1 ; ++$i) {push @$ngram,0}
  if($N>0) {push @$ngram,-1;}

  my $self=
    {
     alphabet=>$alphabet,
     ngram=>$ngram,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $string=$ngramIterator->nextString(); # returns undef if no more
sub nextString
  {
    my ($self)=@_;
    my $alphabet=$self->{alphabet};
    my $ngram=$self->{ngram};
    if(!defined($ngram)) {return undef}
    my $len=@$ngram;
    if($len==0) {undef $self->{ngram}; return ""}
    my $alphaSize=length $alphabet;
    my $i;
    for($i=$len-1 ; $i>=0 ; --$i)
      {
	my $index=++$ngram->[$i];
	if($index<$alphaSize)
	  {return $self->ngramToString()}
	$ngram->[$i]=0;
      }
    return undef;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $self->ngramToString();
sub ngramToString
  {
    my ($self)=@_;
    my $ngram=$self->{ngram};
    my $alphabet=$self->{alphabet};
    my $length=@$ngram;
    my $string;
    for(my $i=0 ; $i<$length ; ++$i)
      {$string.=substr($alphabet,$ngram->[$i],1)}
    return $string;
  }


1;

