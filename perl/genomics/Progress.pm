package Progress;
use strict;

##################################################################
# Attributes:
#   totalWork
#   startTime
# Methods
#   $p=new Progress();
#   $p->start($totalWork);
#   ($timeLeft,$percentDone)=$p->getProgress($workDone);
##################################################################

#-----------------------------------------------------------------
#   new();
sub new
  {
    my ($class)=@_;

    my $self={};
    bless $self,$class;

    return $self;
  }
#-----------------------------------------------------------------
#   start(totalWork);
sub start
  {
    my ($self,$totalWork)=@_;
    
    $self->{totalWork}=$totalWork;
    $self->{startTime}=time();
  }
#-----------------------------------------------------------------
#   ($timeLeft,$percentDone)=getProgress($workDone);
sub getProgress
  {
    my ($self,$workDone)=@_;

    if($workDone==0) {return ("? min",0)}

    my $totalWork=$self->{totalWork};
    my $percentDone=int(100*$workDone/$totalWork*10+0.5)/10;

    my $elapsedSec=time()-$self->{startTime};
    my $secPerUnitWork=$elapsedSec/$workDone;
    my $workRemaining=$totalWork-$workDone;
    my $sec=int($secPerUnitWork*$workRemaining);

    if($sec<60)
      { return ("$sec sec",$percentDone) }
    my $min=int($sec/60*10)/10;
    if($min<60)
      { return ("$min min",$percentDone) }
    my $hours=int($min/60*10)/10;
    if($hours<24)
      { return ("$hours hours",$percentDone) }
    my $days=int($hours/24*10)/10;
    return ("$days days",$percentDone);
  }
#-----------------------------------------------------------------
#-----------------------------------------------------------------


1;
