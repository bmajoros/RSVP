package Fastb;
use strict;
use FastbTrack;

######################################################################
#
# Fastb.pm bmajoros@duke.edu 4/6/2012
#
# Loads a Fastb file into memory and provides access to tracks.
#
# Attributes:
#   trackHash : hash mapping name to FastbTrack
#   trackArray : array of FastbTrack
# Methods:
#   $fastb=new Fastb($filename);
#   $n=$fastb->numTracks();
#   $L=$fastb->getLength();
#   $track=$fastb->getIthTrack($i);
#   $track=$fastb->getTrackByName($id);
#   $fastb->renameTrack($oldName,$newName);
#   $fastb->addTrack($fastbTrack);
#   $fastb->save($filename);
# Private:
#   load($filename);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $fastb=new Fastb($filename);
sub new
{
  my ($class,$filename)=@_;

  my $self=
    {
     trackArray=>[],
     trackHash=>{},
    };
  bless $self,$class;
  $self->load($filename);

  return $self;
}
#---------------------------------------------------------------------
#   $n=$fastb->numTracks();
sub numTracks {
  my ($self)=@_;
  return 0+@{$self->{trackArray}};
}
#---------------------------------------------------------------------
#   $track=$fastb->getIthTrack($i);
sub getIthTrack {
  my ($self,$i)=@_;
  return $self->{trackArray}->[$i];
}
#---------------------------------------------------------------------
#   $track=$fastb->getTrackByName($id);
sub getTrackByName {
  my ($self,$id)=@_;
  return $self->{trackHash}->{$id};
}
#---------------------------------------------------------------------
#   $fastb->addTrack($fastbTrack);
sub addTrack {
  my ($self,$track)=@_;
  push @{$self->{trackArray}},$track;
  $self->{trackHash}->{$track->getID()}=$track;
}
#---------------------------------------------------------------------
#   $fastb->save($filename);
sub save {
  my ($self,$filename)=@_;
  open(OUT,">$filename") || die "can't write to file: $filename\n";
  my $N=$self->numTracks();
  for(my $i=0 ; $i<$N ; ++$i) {
    $self->getIthTrack($i)->save(\*OUT);
  }
  close(OUT);
}
#---------------------------------------------------------------------
#   $fastb->renameTrack($oldName,$newName);
sub renameTrack {
  my ($self,$oldName,$newName)=@_;
  my $hash=$self->{trackHash};
  my $oldTrack=$hash->{$oldName};
  if(!$oldTrack) { return }
  $oldTrack->rename($newName);
  $hash->{$newName}=$oldTrack;
  $hash->{$oldName}=undef;
  delete $hash->{$oldName};
}
#---------------------------------------------------------------------
#   $L=$fastb->getLength();
sub getLength {
  my ($self)=@_;
  if($self->numTracks()==0) { return 0 }
  return $self->getIthTrack(0)->getLength();
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   load($filename);
sub load {
  my ($self,$filename)=@_;
  my @lines;
  open(IN,$filename) || die "Can't open file: $filename\n";
  while(<IN>) {
    chomp;
    if(/\S/) { push @lines,$_ }
  }
  close(IN);
  my $numLines=@lines;
  for(my $i=0 ; $i<$numLines ; ++$i) {
    my $line=$lines[$i];
    #print ">>$line<<\n";
    if($line=~/^\s*([%>])\s*(\S+)(.*)/) {
      my ($op,$id,$rest)=($1,$2,$3);
      if($op eq ">") {
	my $seq;
	for(++$i ; $i<$numLines ; ++$i) {
	  $line=$lines[$i];
	  if($line=~/^\s*([%>])\s*(\S+)(.*)/) { --$i; last }
	  $seq.=$line;
	}
	my $track=new FastbTrack("discrete",$id,$seq,$rest);
	$self->addTrack($track);
      }
      else {
	my $data=[];
	for(++$i ; $i<$numLines ; ++$i) {
	  $line=$lines[$i];
	  if($line=~/^\s*([%>])\s*(\S+)(.*)/) { --$i; last }
	  $line=~/(\S+)/;
	  push @$data,$1;
	}
	my $track=new FastbTrack("continuous",$id,$data,$rest);
	$self->addTrack($track);
      }
    }
    else { die "can't parse line: $line\n" }
  }
}
#---------------------------------------------------------------------

1;

