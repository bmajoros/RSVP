package BandedSmithWaterman;
use strict;
use NeedlemanAlignment;

######################################################################
#
# BandedSmithWaterman.pm bmajoros@tigr.org 5/17/2004
#
# Runs the "banded-smith-waterman" program to do alignment.  Returns 
# an object of type NeedlemanAlignment.
#
# Attributes:
#   executableDir    : i.e., "/home/bmajoros/GF/twain/alignment"
#   matrixFileNoPath : i.e., "blosum62"
# Methods:
#   $aligner=new BandedSmithWaterman($executableDir,$matrixFileNoPath);
#   $alignment=$aligner->align($fasta1,$fasta2,"DNA|PROTEIN",
#                              $gapOpenPenalty,$gapExtendPenalty,
#                              $bandWidth);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
# $aligner=new BandedSmithWaterman($executableDir,$matrixFileNoPath);
sub new
{
  my ($class,$executableDir,$matrixFileNoPath)=@_;

  if($executableDir=~/\/$/) {chop $executableDir}
  my $self=
    {
     executableDir=>$executableDir,
     matrixFileNoPath=>$matrixFileNoPath,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
# $alignment=$aligner->align($fasta1,$fasta2,"DNA",$gapOpen,$gapExtend,
#                            $bandWidth);
sub align
  {
    my ($self,$fasta1,$fasta2,$type,$gapOpenPenalty,$gapExtendPenalty,
	$bandWidth)=@_;
    my $executableDir=$self->{executableDir};
    my $matrixFileNoPath=$self->{matrixFileNoPath};
    my $exec="$executableDir/banded-smith-waterman";
    my $matrix="$executableDir/matrices/$matrixFileNoPath";
    my $command=
      "$exec $matrix $gapOpenPenalty $gapExtendPenalty $fasta1 $fasta2 "
	."$type $bandWidth";
    my ($identity,$matches,$mismatches,$insertions,$length,$score);
    my (@alignment);
    open(PIPE,"$command|") || die "Can't pipe in from: \"$command\"";
    while(<PIPE>)
      {
       if(/Percent identity:\s*(\d+)%,\s*matches=(\d+),\s*mismatches=(\d+)/)
	  {
	    ($identity,$matches,$mismatches)=($1,$2,$3);
	  }
	elsif(/insertions=(\d+),\s*alignment length=(\d+),\s*score=(\d+)/)
	  {
	    ($insertions,$length,$score)=($1,$2,$3);
	  }
	elsif(/Query: (.*)/)
	  {
	    $alignment[0].=$1;
	    $_=<PIPE>;
	    $_=~/       (.*)/;
	    $alignment[1].=$1;
	    $_=<PIPE>;
	    $_=~/Sbjct: (.*)/;
	    $alignment[2].=$1;
	  }
      }
    close(PIPE);
    return new NeedlemanAlignment(\@alignment,$fasta1,$fasta2,
				  $identity,$matches,$mismatches,
				  $insertions,$length,$score);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

