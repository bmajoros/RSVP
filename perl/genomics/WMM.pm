package WMM;
use strict;

######################################################################
#
# WMM.pm bmajorostigr.org 8/12/2005
#
# Computes log scores of a weight matrix at a given position in a
# DNA sequence.  Loads itself from GeneZilla model files.
#
# Attributes:
#   matrix = array of hashes, each mapping chars to log(P) values
#   length = number of positions in the WMM
#   cutoff = threshold value for significant matches of the WMM
# Methods:
#   $wmm=new WMM($filename);
#   $score=$wmm->score($seq,$pos);
#   $length=$wmm->getLength();
#   $cutoff=$wmm->getCutoff();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $wmm=new WMM($filename);
sub new
{
  my ($class,$filename)=@_;

  my $matrix=[];
  my $self=
    {
     cutoff=>0,
     matrix=>$matrix,
    };
  bless $self,$class;

  $self->load($filename);

  return $self;
}
#---------------------------------------------------------------------
#   $score=$wmm->score($seq,$pos);
sub score
  {
    my ($self,$seq,$pos)=@_;
    my $length=$self->getLength();
    my $matrix=$self->{matrix};
    my $score=0;
    for(my $i=0 ; $i<$length; ++$i)
      {
	my $s=substr($seq,$pos+$i,1);
	my $logP=$matrix->[$i]->{$s};
	$score+=$logP;
      }
    return $score;
  }
#---------------------------------------------------------------------
#   $length=$wmm->getLength();
sub getLength
  {
    my ($self)=@_;
    return $self->{length};
  }
#---------------------------------------------------------------------
#   $cutoff=$wmm->getCutoff();
sub getCutoff
  {
    my ($self)=@_;
    return $self->{cutoff};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
sub load
  {
    my ($self,$filename)=@_;
    my $msg="Syntax error in model file ($filename)\n";
    open(IN,$filename) || die "Can't open file: $filename\n";
    my $line=<IN>;
    die $msg unless $line=~/^WMM\s*$/;
    $line=<IN>; # signa type
    $line=<IN>;
    $line=~/^(\S+)\s+(\d+)\s+\d+\s*$/ || die $msg;
    my ($cutoff,$length)=($1,$2);
    $self->{cutoff}=$cutoff;
    $self->{length}=$length;
    <IN>;
    <IN>; # strand
    my $matrix=$self->{matrix};
    my @symbols=('A','C','G','N','T');
    for(my $i=0 ; $i<$length ; ++$i)
      {
	$line=<IN>;
	my @fields=split/\s+/,$line;
	for(my $j=0 ; $j<5 ; ++$j)
	  {
	    my $s=$symbols[$j];
	    my $logP=$fields[$j];
	    if($logP eq "-inf") {$logP=log(0.000000000000001)}
	    $matrix->[$i]->{$s}=$logP;
	  }
      }
    close(IN);
  }


1;

