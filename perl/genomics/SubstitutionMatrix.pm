package SubstitutionMatrix;
use strict;
use Alphabet;

######################################################################
#
# SubstitutionMatrix.pm bmajorostigr.org 10/9/2003
#
# 
# 
#
# Attributes:
#
# Methods:
#   $matrix=new SubstitutionMatrix($filename);
#   $alphabet=$matrix->getAlphabet();
#   my $value=$matrix->lookup($residue1,$residue2);
#   my $newResidue=$matrix->chooseLikelySubstitution($oldResidue);
# Private methods:
#   $self->load($filename);
#   $self->computeProbTable();
#   
######################################################################

my $PSEUDO_COUNTS=0.05;

#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $substitutionMatrix=new SubstitutionMatrix($filename);
sub new
{
  my ($class,$filename)=@_;
  
  my $self={};
  bless $self,$class;

  $self->load($filename);
  $self->computeProbTable();

  return $self;
}
#---------------------------------------------------------------------
#   my $value=$matrix->lookup($residue1,$residue2);
sub lookup
  {
    my ($self,$residue1,$residue2)=@_;
    my $matrix=$self->{matrix};
    my $alphabet=$self->{alphabet};
    my $index1=$alphabet->symbolToIndex($residue1);
    my $index2=$alphabet->symbolToIndex($residue2);
    return $matrix->[$index1][$index2];
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
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
    open(IN,$filename) || die "Can't open $filename in SubstitutionMatrix::load()";
    my $matrix=[];
    my $alphabet;
    while(<IN>)
      {
	next if(/\#/);
	if(/^\s+(.*)/)
	  {
	    my $header=$1;
	    $header=~s/\s+//g;
	    $alphabet=new Alphabet($header);
	  }
	else
	  {
	    my @fields=split/\s+/,$_;
	    my $rowHeader=shift @fields;
	    my $index=$alphabet->symbolToIndex($rowHeader);
	    $matrix->[$index]=\@fields;
	  }
      }
    close(IN);
    $self->{matrix}=$matrix;
    $self->{alphabet}=$alphabet;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#   my $newResidue=$matrix->chooseLikelySubstitution($oldResidue);
sub chooseLikelySubstitution
  {
    my ($self,$oldResidue)=@_;
    my $probTable=$self->{probTable};
    my $alphabet=$self->{alphabet};
    my $nAlpha=$alphabet->size();
    my $p=rand(1);
    my $i;
    for($i=0 ; $i<$nAlpha ; ++$i)
      {
	my $upperBoundary=$probTable->[$oldResidue][$i];
	last unless $p>$upperBoundary;
      }
    if($i==$nAlpha) {--$i}
    my $oldIndex=$alphabet->symbolToIndex($oldResidue);
    if($i==$oldIndex) {return $self->chooseLikelySubstitution($oldResidue)}
    return $alphabet->indexToSymbol($i);
  }
#---------------------------------------------------------------------
#   $self->computeProbTable();
sub computeProbTable
  {
    my ($self)=@_;
    my $matrix=$self->{matrix};
    my $alphabet=$self->{alphabet};
    my $numSymbols=$alphabet->size();
    my $probTable=[];
    for(my $x=0 ; $x<$numSymbols ; ++$x)
      {
	my $min=0;
	for(my $y=0 ; $y<$numSymbols ; ++$y)
	  {
	    my $value=$matrix->[$x][$y];
	    $min=($value<$min ? $value : $min);
	  }
	my $offset=($min<0 ? -$min : 0);
	my $sum=0;
	for(my $y=0 ; $y<$numSymbols ; ++$y)
	  {
	    my $value=$matrix->[$x][$y]+$offset+$PSEUDO_COUNTS;
	    $sum+=$value;
	  }
	my $accum=0;
	for(my $y=0 ; $y<$numSymbols ; ++$y)
	  {
	    my $P=($matrix->[$x][$y]+$offset+$PSEUDO_COUNTS)/$sum;
	    $accum+=$P;
	    $probTable->[$x][$y]=$accum;
	  }
      }
    $self->{probTable}=$probTable;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------


1;

