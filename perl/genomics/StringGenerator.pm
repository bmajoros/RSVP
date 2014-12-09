package StringGenerator;
use strict;
use Carp;

######################################################################
#
# StringGenerator.pm bmajoros
#
# Generates random strings of specified length, over a given alphabet,
# and according to a given distribution.  If you specify just the
# number of symbols, it will use a uniform distribution over an
# alphabet of arbitrary symbols.
#
# Attributes:
#   alphabet : array of symbols
#   distribution : hash table mapping symbols to probabilities
#   rouletteWheel
# Methods:
#   $g=new StringGenerator({numSymbols=>$n});
#   $g=new StringGenerator({distribution=>{a=>.3,b=>.4,c=>.2}});
#   $stringReference=$g->generateString($length);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $g=new StringGenerator({numSymbols=>$n});
#   $g=new StringGenerator({distribution=>{a=>.3,b=>.4,c=>.2}});
sub new
{
  my ($class,$parms)=@_;
  
  my $self=
    {
     alphabet=>[],
     distribution=>{},
    };
  bless $self,$class;

  if(defined($parms->{numSymbols})) 
    { 
      $self->initUniform($parms->{numSymbols}) 
    }
  elsif(defined($parms->{distribution})) 
    {
      my $distribution=$parms->{distribution};
      my @alphabet=keys %$distribution;
      $self->{alphabet}=\@alphabet;
      $self->{distribution}=$distribution;
    }
  else
    {
      confess "StringGenerator::new() called without proper arguments";
    }

  $self->initWheel();

  return $self;
}
#---------------------------------------------------------------------
#   $string=$g->generateString($length);
sub generateString
  {
    my ($self,$length)=@_;

    my $string;
    for(my $i=0 ; $i<$length ; ++$i)
      {
	my $p=rand();
	my $symbol=$self->findSymbol($p);
	$string.=$symbol;
      }
    return \$string;
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $self->initUniform($numSymbols);
sub initUniform
  {
    my ($self,$numSymbols)=@_;

    my $a=unpack("c","a");
    my @alphabet;
    my %distribution;
    my $p=1/$numSymbols;
    for(my $i=0 ; $i<$numSymbols ; ++$i)
      {
	my $symbol=pack("c",$a+$i);
	push @alphabet,$symbol;
	$distribution{$symbol}=$p;
      }

    $self->{alphabet}=\@alphabet;
    $self->{distribution}=\%distribution;
  }
#---------------------------------------------------------------------
sub initWheel
  {
    my ($self)=@_;

    my @wheel;
    my $alphabet=$self->{alphabet};
    my $numSymbols=@$alphabet;
    my $distribution=$self->{distribution};
    my $p=0;
    for(my $i=0 ; $i<$numSymbols ; ++$i)
      {
	my $symbol=$alphabet->[$i];
	my $slice=$distribution->{$symbol};

	$wheel[$i][0]=$p;
	$wheel[$i][1]=$symbol;
	$p+=$slice;
      }
    $self->{rouletteWheel}=\@wheel;
  }
#---------------------------------------------------------------------
sub findSymbol
  {
    my ($self,$p)=@_;
    
    my $wheel=$self->{rouletteWheel};
    my $n=@$wheel;
    my $begin=0;
    my $end=$n;
    while($begin<$end-1)
      {
	my $mid=int(($begin+$end)/2);
	my $midValue=$wheel->[$mid][0];
	if($p>$midValue) { $begin=$mid }
	else { $end=$mid }
      }
    return $wheel->[$begin][1];
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------



1;

