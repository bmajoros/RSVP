package Interval;
use strict;

######################################################################
#
# Interval.pm bmajoros@duke.edu 4/12/2006
#
# Represents a range [b,e) where b is included in the range and e is
# not.  Any interval of the form [x,x) is considered empty.  The length
# of an interval [b,e) is defined as e-b.  
#
# Attributes:
#   begin : index of the first element in the interval
#   end   : index of the element AFTER the last element in the interval
# Methods:
#   $interval=new Interval($begin,$end);
#   $begin=$interval->getBegin();
#   $end=$interval->getEnd();
#   $bool=$interval->overlaps($other);
#   $bool=$interval->contains($index);
#   $intersection=$interval->intersect($other);
#   $union=$interval->union($other);
#   $diff=$interval->minus($other);  # returns a list of intervals!
#   $length=$interval->getLength();
#   $bool=$interval->equals($other);
#   $other=$interval->clone();
#   $bool=$interval->isEmpty();
#   $d=$interval->relativeDistanceFromBegin($pos);
#   $d=$interval->relativeDistanceFromEnd($pos);
#   $interval->printOn(\*STDOUT);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $interval=new Interval($begin,$end);
sub new
{
  my ($class,$begin,$end)=@_;
  
  my $self=
  {
    begin=>$begin,
    end=>$end,
  };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $interval->printOn(\*STDOUT);
sub printOn {
    my ($self,$fileHandle)=@_;
    my $begin=$self->getBegin();
    my $end=$self->getEnd();
    print $fileHandle "$begin,$end";
}
#---------------------------------------------------------------------
#   $begin=$interval->getBegin();
sub getBegin
{
    my ($self)=@_;
    return $self->{begin};
}
#---------------------------------------------------------------------
#   $end=$interval->getEnd();
sub getEnd
{
    my ($self)=@_;
    return $self->{end};
}
#---------------------------------------------------------------------
#   $bool=$interval->overlaps($other);
sub overlaps
{
    my ($self,$other)=@_;
    my $begin=$self->{begin};
    my $end=$self->{end};
    my $otherBegin=$other->{begin};
    my $otherEnd=$other->{end};
    return $begin<$otherEnd && $otherBegin<$end;
}
#---------------------------------------------------------------------
#   $intersection=$interval->intersect($other);
sub intersect
{
    my ($self,$other)=@_;
    if($self->isEmpty()) {return $self->clone()}
    if($other->isEmpty()) {return $other->clone()}
    my $begin=$self->{begin};
    my $end=$self->{end};
    my $otherBegin=$other->{begin};
    my $otherEnd=$other->{end};
    if($otherBegin>$begin) {$begin=$otherBegin}
    if($otherEnd<$end) {$end=$otherEnd}
    if($end<$begin) {$end=$begin}
    return new Interval($begin,$end);
}
#---------------------------------------------------------------------
#   $union=$interval->union($other);
sub union
{
    my ($self,$other)=@_;
    if($self->isEmpty()) {return $other->clone()}
    if($other->isEmpty()) {return $self->clone()}
    my $begin=$self->{begin};
    my $end=$self->{end};
    my $otherBegin=$other->{begin};
    my $otherEnd=$other->{end};
    if($otherBegin<$begin) {$begin=$otherBegin}
    if($otherEnd>$end) {$end=$otherEnd}
    return new Interval($begin,$end);
}
#---------------------------------------------------------------------
#   $length=$interval->getLength();
sub getLength
{
    my ($self)=@_;
    my $length=$self->{end}-$self->{begin};
    if($length<0) {$length=0}
    return $length;
}
#---------------------------------------------------------------------
#   $bool=$interval->equals($other);
sub equals
{
    my ($self,$other)=@_;
    return 
      $self->isEmpty() && $other->isEmpty()
        ||
      $self->{begin}==$other->{begin} &&
        $self->{end}==$other->{end};
}
#---------------------------------------------------------------------
#   $other=$interval->clone();
sub clone
{
    my ($self)=@_;
    my $other=new Interval($self->{begin},$self->{end});
    return $other;
}
#---------------------------------------------------------------------
#   $diff=$interval->minus($other);
sub minus
{
    my ($self,$other)=@_;
    if(!$self->overlaps($other) || $other->isEmpty() || $self->isEmpty())
      {return [$self->clone()]}
    my $begin=$self->{begin};
    my $end=$self->{end};
    my $otherBegin=$other->{begin};
    my $otherEnd=$other->{end};
    if($otherBegin<=$begin)
    {
        if($otherEnd>=$end)
        {
            return [];
        }
        else
        {
            return [new Interval($otherEnd,$end)];
        }
    }
    else
    {
        if($otherEnd>=$end)
        {
            return [new Interval($begin,$otherBegin)];
        }
        else
        {
            my $firstInterval=new Interval($begin,$otherBegin);
            my $secondInterval=new Interval($otherEnd,$end);
            return [$firstInterval,$secondInterval];
        }
    }
}
#---------------------------------------------------------------------
#   $bool=$interval->isEmpty();
sub isEmpty
{
    my ($self)=@_;
    return $self->{begin}>=$self->{end};
}
#---------------------------------------------------------------------
#   $bool=$interval->contains($index);
sub contains
{
    my ($self,$index)=@_;
    return $index>=$self->{begin} && $index<$self->{end};
}
#---------------------------------------------------------------------
#   $d=$interval->relativeDistanceFromBegin($pos);
sub relativeDistanceFromBegin
{
    my ($self,$pos)=@_;
    if(!$self->contains($pos)) {die}
    my $L=$self->getLength();
    return ($pos-$self->{begin})/$L;
}
#---------------------------------------------------------------------
#   $d=$interval->relativeDistanceFromEnd($pos);
sub relativeDistanceFromEnd
{
    my ($self,$pos)=@_;
    if(!$self->contains($pos)) {die}
    my $L=$self->getLength();
    my $revPos=$self->{end}-$pos-1;
    return $revPos/$L;
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

