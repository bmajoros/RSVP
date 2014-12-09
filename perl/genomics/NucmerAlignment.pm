package NucmerAlignment;
use strict;
use NucmerHit;

######################################################################
#
# NucmerAlignment.pm bmajoros@tigr.org 4/26/2005
#
#   A collection of NucmerHits.
#
# Attributes:
#   hits : array of NucmerHit objects
#
# Methods:
#   $alignment=new NucmerAlignment($deltaFilename);
#   $n=$alignment->getNumHits();
#   $hit=$alignment->getIthHit($i);
# Private methods:
#   $alignment->load($filename);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class,$filename)=@_;

  my $hits=[];
  my $self=
    {
     hits=>$hits
    };
  bless $self,$class;

  $self->load($filename);

  return $self;
}
#---------------------------------------------------------------------
#   $n=$alignment->getNumHits();
sub getNumHits
  {
    my ($self)=@_;
    my $hits=$self->{hits};
    my $numHits=@$hits;
    return $numHits;
  }
#---------------------------------------------------------------------
#   $hit=$alignment->getIthHit($i);
sub getIthHit
  {
    my ($self,$i)=@_;
    my $hits=$self->{hits};
    return $hits->[$i];
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $alignment->load($filename);
sub load
  {
    my ($self,$filename)=@_;
    my $hits=$self->{hits};
    open(IN,$filename) || die "Can't open file: $filename\n";
    <IN>;<IN>;<IN>;
    while(1)
      {
	my $hit=new NucmerHit(\*IN);
	last unless defined $hit;
	push @$hits,$hit;
      }
    close(IN);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


1;

