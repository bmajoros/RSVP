package AppalFile;
use strict;

######################################################################
#
# AppalFile.pm bmajorostigr.org 4/26/2005
#
# Represents an alignment produced by the APPAL program (part of TWAIN)
# and stored in a *.appal file.
#
# Attributes:
#   cells
#
# Methods:
#   $appalFile=new AppalFile($filename);
#   $numCells=$appalFile->getNumCells();
#   $cell=$appalFile->getIthCell($i);
#
# Private Methods:
#   $self->load($filename);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $appalFile=new AppalFile($filename);
sub new
{
  my ($class,$filename)=@_;

  my $cells=[];
  my $self={cells=>$cells};
  bless $self,$class;

  $self->load($filename);

  return $self;
}
#---------------------------------------------------------------------
#   $numCells=$appalFile->getNumCells();
sub getNumCells
  {
    my ($self)=@_;
    my $cells=$self->{cells};
    return 0+@$cells;
  }
#---------------------------------------------------------------------
#   $cell=$appalFile->getIthCell($i);
sub getIthCell
  {
    my ($self,$i)=@_;
    return $self->{cells}->[$i];
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $self->load($filename);
sub load
  {
    my ($self,$filename)=@_;
    my $cells=$self->{cells};
    open(IN,$filename) || die "Can't open file: $filename\n";
    my $header=<IN>;
    my ($qBegin,$qEnd,$rBegin,$rEnd)=split/\s+/,$header;
    --$qBegin;
    --$qEnd;
    my $qPos=$qBegin;
    my $rPos=$rBegin;
    push @$cells,[$qPos,$rPos];
    while(<IN>)
      {
	if(/(\S)\s+(\d+)/)
	  {
	    my ($op,$len)=($1,$2);
	    if($op eq "i")
	      {
		for(my $i=0 ; $i<$len ; ++$i)
		  {push @$cells,[$qPos,++$rPos]}
	      }
	    elsif($op eq "d")
	      {
		for(my $i=0 ; $i<$len ; ++$i)
		  {push @$cells,[++$qPos,$rPos]}
	      }
	    elsif($op eq "m")
	      {
		for(my $i=0 ; $i<$len ; ++$i)
		  {push @$cells,[++$qPos,++$rPos]}
	      }
	    elsif($op eq "e")
	      { last }
	    else {die "bad operator: $op"}
	  }
      }
    close(IN);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


1;

