package NucmerHit;
use strict;
use FileHandle;

######################################################################
#
# NucmerHit.pm bmajoros@tigr.org 4/26/2005
#
# 
# 
#
# Attributes:
#  rBegin
#  rEnd
#  qBegin
#  qEnd
#  cells
#  deltas
#
# Methods:
#   $nucmerHit=new NucmerHit(rBegin,rEnd,qBegin,qEnd,cells);
#   $cells=$alignment->getCells();
#   $rBegin=$alignment->getRefBegin();
#   $rEnd=$alignment->getRefEnd();
#   $qBegin=$alignment->getQueryBegin();
#   $qEnd=$alignment->getQueryEnd();
#
# Private methods:
#   $hit->load($filehandle);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class,$filehandle)=@_;

  my $cells=[];
  my $self=
    {
     cells=>$cells
    };
  bless $self,$class;

  return $self->load($filehandle);
}
#---------------------------------------------------------------------
#   $cells=$alignment->getCells();
sub getCells
  {
    my ($self)=@_;
    return $self->{cells};
  }
#---------------------------------------------------------------------
#   $rBegin=$alignment->getRefBegin();
sub getRefBegin
  {
    my ($self)=@_;
    return $self->{rBegin};
  }
#---------------------------------------------------------------------
#   $rEnd=$alignment->getRefEnd();
sub getRefEnd
  {
    my ($self)=@_;
    return $self->{rEnd};
  }
#---------------------------------------------------------------------
#   $qBegin=$alignment->getQueryBegin();
sub getQueryBegin
  {
    my ($self)=@_;
    return $self->{qBegin};
  }
#---------------------------------------------------------------------
#   $qEnd=$alignment->getQueryEnd();
sub getQueryEnd
  {
    my ($self)=@_;
    return $self->{qEnd};
  }
#---------------------------------------------------------------------
sub load
  {
    my ($self,$filehandle)=@_;

    # First, load the header and the deltas:
    my $header=<$filehandle>;
    return undef unless defined($header);
    my @fields=split/\s+/,$header;
    return undef unless @fields==7;
    my ($queryBegin,$queryEnd,$refBegin,$refEnd,$numIdent,
	$numSim,$numStops)=@fields;
    --$refBegin; --$queryBegin;
    my $refStrand=$refBegin<$refEnd ? '+' : '-';
    my $queryStrand=$queryBegin<$queryEnd ? '+' : '-';
    if($refStrand eq '-') 
      {($refBegin,$refEnd)=($refEnd,$refBegin)}
    if($queryStrand eq '-') 
      {($queryBegin,$queryEnd)=($queryEnd,$queryBegin)}
    $self->{qBegin}=$queryBegin;
    $self->{qEnd}=$queryEnd;
    $self->{qStrand}=$queryStrand;
    $self->{rBegin}=$refBegin;
    $self->{rEnd}=$refEnd;
    $self->{rStrand}=$refStrand;
    $self->{numIdent}=$numIdent;
    my @deltas;
    while(1)
      {
	my $delta=0+<$filehandle>;
	last unless $delta;
	push @deltas,$delta;
      }
    $self->{deltas}=\@deltas;

    # Now convert the deltas into an array of cells:
    my $cells=$self->{cells};
    my $numDeltas=@deltas;
    my $rPos=$refBegin;
    my $qPos=$queryBegin;
    push @$cells,[$rPos,$qPos];
    for(my $i=0 ; $i<$numDeltas ; ++$i)
      {
	my $delta=$deltas[$i];
	my $absDelta=abs($delta)-1;
	for(my $d=0 ; $d<$absDelta ; ++$d)
	  {
	    ++$rPos;
	    ++$qPos;
	    push @$cells,[$rPos,$qPos];
	  }
	if($delta>0) {++$qPos}
	else {++$rPos}
	push @$cells,[$rPos,$qPos];
      }
    while($qPos<$queryEnd-1)
      {
	++$qPos;
	++$rPos;
	push @$cells,[$rPos,$qPos];
      }

    return $self;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

1;

