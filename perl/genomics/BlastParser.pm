package BlastParser;
use strict;
$|=1;

######################################################################
#
# BlastParser.pm bmajorostigr.org 12/12/2003
#
# Parses a blast report and returns a list of pairs, where the first
# element of the pair is the query ID and the second element is an
# object: {id=>$id,defline=>$def,identity=>$ident,expect=>$exp,
# alignLength=>$alignLength}.  The pairs are sorted by E-value.
#
# Attributes:
#
# Methods:
#   $blastParser=new BlastParser();
#   $hits=$blastParser->parse($filename,$maxExpect);
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self={};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $hits=$blastParser->parse($filename,$maxExpect);
sub parse
  {
    my ($self,$filename,$maxExpect)=@_;
    my $hash={};
    open(IN,$filename) || die "Can't open file: $filename\n";
    my ($queryName,$queryLength,$subjectLength,$defline,$expect,$identity,
	$subjectId,$alignLength,$queryDefline);
    while(<IN>)
      {
	if(/Query=(.*)/)
	  {
	    $queryName=$1;
	    while(<IN>)
	      {
		if(/\((\d+)\s+letters/)
		  {
		    $queryLength=$1;
		    chomp;
		    $queryName.=" $_";
		    last;
		  }
	      }
	    if($queryName=~/^\s*([^;\s]+)(.*)/) 
	      {
		$queryName=$1;
		$queryDefline="$1$2";
	      }
	  }
	elsif(/^>(.*)/)
	  {
	    $defline=$1;
	    if($defline=~/^([^;\s]+)\s*(.*)/) 
	      {$subjectId=$1;$defline=$2} 
	    else {undef $subjectId}
	    while(<IN>)
	      {
		if(/Length\s*=\s*(\d+)/)
		  {
		    $subjectLength=$1;
		    last;
		  }
		chomp;
		$defline.=" $_";
	      }
	  }
	elsif(/Score\s*=.*Expect\s*=\s*([^,;\s]+)/)
	  {
	    $expect=$1;
	  }
	elsif(/Identities\s*=\s*(\d+)\/(\d+)/)
	  {
	    $alignLength=$2;
	    my $identities=$1;
	    $identity=$identities/$alignLength;
	    if($expect<=$maxExpect)
	      {
		if(!defined($hash->{$queryName}))
		  {$hash->{$queryName}=[]}
		my $array=$hash->{$queryName};
		push @$array,
		  {
		   id=>$subjectId,
		   defline=>$defline,
		   identity=>$identity,
		   expect=>$expect,
		   alignLength=>$alignLength,
		   queryLength=>$queryLength,
		   identities=>$identities
		  };
	      }
	  }
      }
    close(IN);

    my $hits=[];
    my @keys=keys %$hash;
    my $n=@keys;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $key=$keys[$i];
	my $matches=$hash->{$key};
	foreach my $match(@$matches)
	  {
	    push @$hits,[$key,$match];
	  }
      }
    @$hits=sort {$a->[1]->{expect} <=> $b->[1]->{expect}} @$hits;
    return $hits;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

