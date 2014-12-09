package FastbTrack;
use strict;
use FastaWriter;

######################################################################
#
# FastbTrack.pm bmajoros@duke.edu 4/6/2012
#
# A track in a Fastb file.
#
# Attributes:
#   type : "discrete" or "continuous"
#   id : name of track
#   data : string (for discrete) or array of float (for continuous)
#   deflineExtra : extra info for defline
# Methods:
#   $track=new FastbTrack($type,$id,$data,$deflineExtra); # $type="discrete" or "continuous"
#   $type=$track->getType();
#   $data=$track->getData();
#   $id=$track->getID();
#   $L=$track->getLength();
#   $track->rename($newID);
#   $bool=$track->isDiscrete();
#   $bool=$track->isContinuous();
#   $track->save(FILEHANDLE);
#   $array=$track->getNonzeroRegions(); # returns array of [begin,end]
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $track=new FastbTrack($type,$id,$data,$deflineExtra); # $type="discrete" or "continuous"
sub new
{
  my ($class,$type,$id,$data,$deflineExtra)=@_;
  my $self=
    {
     type=>$type,
     id=>$id,
     data=>$data,
     deflineExtra=>$deflineExtra
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $type=$track->getType();
sub getType  {
  my ($self)=@_;
  return $self->{type};
}
#---------------------------------------------------------------------
#   $data=$track->getData();
sub getData  {
  my ($self)=@_;
  return $self->{data};
}
#---------------------------------------------------------------------
#   $bool=$track->isDiscrete();
sub isDiscrete  {
  my ($self)=@_;
  return $self->{type} eq "discrete";
}
#---------------------------------------------------------------------
#   $bool=$track->isContinuous();
sub isContinuous {
  my ($self)=@_;
  return $self->{type} eq "continuous";
}
#---------------------------------------------------------------------
#   $id=$track->getID();
sub getID {
  my ($self)=@_;
  return $self->{id};
}
#---------------------------------------------------------------------
#   $track->save(FILEHANDLE);
sub save {
  my ($self,$fh)=@_;
  my $writer=new FastaWriter;
  my $id=$self->{id};
  my $data=$self->{data};
  my $deflineExtra=$self->{deflineExtra};
  if($self->isDiscrete()) { $writer->addToFasta(">$id $deflineExtra",
						$data,$fh) }
  else {
    print $fh "%$id $deflineExtra\n";
    my $n=@$data;
    for(my $i=0 ; $i<$n ; ++$i) {
      my $line=$data->[$i];
      print $fh "$line\n";
    }
  }
}
#---------------------------------------------------------------------
#   $array=$track->getNonzeroRegions(); # returns array of [begin,end]
sub getNonzeroRegions {
  my ($self)=@_;
  my $data=$self->{data};
  if($self->isDiscrete()) { die "track is not continuous" }
  my $L=@$data;
  my $intervals=[];
  my $begin;
  for(my $i=0 ; $i<$L ; ++$i) {
    my $x=$data->[$i];
    if($x>0 && ($i==0 || $data->[$i-1]==0)) { $begin=$i }
    elsif($x==0 && $i>0 && $data->[$i-1]>0) {
      push @$intervals,[$begin,$i];
    }
  }
  if($L>0 && $data->[$L-1]>0) {
    push @$intervals,[$begin,$L];
  }
  return $intervals;
}
#---------------------------------------------------------------------
#   $track->rename($newID);
sub rename {
  my ($self,$newID)=@_;
  $self->{id}=$newID;
}
#---------------------------------------------------------------------
#   $L=$track->getLength();
sub getLength {
  my ($self)=@_;
  my $data=$self->getData();
  if($self->isDiscrete()) { return length $data }
  else { return 0+@$data }
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------







#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

