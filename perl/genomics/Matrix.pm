package Matrix;
use strict;

######################################################################
#
# Matrix.pm bmajoros@duke.edu 11/17/2005
#
#
# Attributes:
#   array : 2D array that stores the elements
#   numRows
#   numColumns
# Methods:
#   $matrix=new Matrix($m,$n);
#   $m=$matrix->numRows();
#   $n=$matrix->numColumns();
#   $v=$matrix->get($i,$j);
#   $matrix->set($i,$j,$value);
#   $matrix->save($filehandle);
#   $matrix->load($filehandle);
#   $row=$matrix->getRow($i);
#   $col=$matrix->getColumn($j);
#   $s=$matrix->submatrix($fromRow,$toRow,$fromCol,$toCol); #inclusive!
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $matrix=new Matrix($m,$n);
sub new
{
  my ($class,$m,$n)=@_;

  $m+=0;
  $n+=0;
  my $array=[];
  for(my $i=0 ; $i<$m ; ++$i)
    {
      for(my $j=0 ; $j<$n ; ++$j)
	{
	  $array->[$i]->[$j]=0;
	}
    }

  my $self=
    {
     array=>$array,
     numRows=>$m,
     numColumns=>$n
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $m=$matrix->numRows();
sub numRows
  {
    my ($self)=@_;
    return $self->{numRows};
  }
#---------------------------------------------------------------------
#   $n=$matrix->numColumns();
sub numColumns
  {
    my ($self)=@_;
    return $self->{numColumns};
  }
#---------------------------------------------------------------------
#   $v=$matrix->get($i,$j);
sub get
  {
    my ($self,$i,$j)=@_;
    return $self->{array}->[$i]->[$j];
  }
#---------------------------------------------------------------------
#   $matrix->set($i,$j,$value);
sub set
  {
    my ($self,$i,$j,$value)=@_;
    $self->{array}->[$i]->[$j]=$value;
  }
#---------------------------------------------------------------------
#   $matrix->save($filehandle);
sub save
  {
    my ($self,$fh)=@_;
    my $m=$self->{numRows};
    my $n=$self->{numColumns};
    print $fh "$m\t$n\n";
    my $array=$self->{array};
    for(my $i=0 ; $i<$m ; ++$i)
      {
	my $row=$array->[$i];
	for(my $j=0 ; $j<$n ; ++$j)
	  {
	    my $x=$row->[$j];
	    print $fh "$x\t";
	  }
	print $fh "\n";
      }
  }
#---------------------------------------------------------------------
#   $matrix->load($filehandle);
sub load
  {
    my ($self,$fh)=@_;
    my $array=[];
    $_=<$fh>;
    my @fields=split/\s+/,$_;
    my ($m,$n)=@fields;
    for(my $i=0 ; $i<$m ; ++$i)
      {
	$_=<$fh>;
	@fields=split/\s+/,$_;
	for(my $j=0 ; $j<$n ; ++$j)
	  {
	    $array->[$i]->[$j]=0+$fields[$j];
	  }
      }
    $self->{numRows}=$m;
    $self->{numColumns}=$n;
    $self->{array}=$array;
  }
#---------------------------------------------------------------------
#   $row=$matrix->getRow($i);
sub getRow
  {
    my ($self,$i)=@_;
    my $n=$self->{numColumns};
    my $matrix=new Matrix(1,$n);
    for(my $j=0 ; $j<$n ; ++$j)
      {
	my $value=$self->get($i,$j);
	$matrix->set($i,$j,$value);
      }
    return $matrix;
  }
#---------------------------------------------------------------------
#   $col=$matrix->getColumn($j);
sub getColumn
  {
    my ($self,$j)=@_;
    my $m=$self->{numRows};
    my $matrix=new Matrix($m,1);
    for(my $i=0 ; $i<$m ; ++$i)
      {
	my $value=$self->get($i,$j);
	$matrix->set($i,$j,$value);
      }
    return $matrix;
  }
#---------------------------------------------------------------------
#   $s=$matrix->submatrix($fromRow,$toRow,$fromCol,$toCol);
sub submatrix
  {
    my ($self,$fromRow,$toRow,$fromCol,$toCol)=@_;
    my $m=$toRow-$fromRow+1;
    my $n=$toCol-$fromCol+1;
    my $matrix=new Matrix($m,$n);
    for(my $i=$fromRow ; $i<=$toRow ; ++$i)
      {
	for(my $j=$fromCol ; $j<=$toCol ; ++$j)
	  {
	    my $x=$self->get($i,$j);
	    $matrix->set($i-$fromRow,$j-$fromCol,$x);
	  }
      }
    return $matrix;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------



#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

