package SmithWaterman;
use strict;
use NeedlemanAlignment;
use TempFilename;
use FastaWriter;

######################################################################
#
# SmithWaterman.pm bmajoros@tigr.org 5/17/2004
#
# Runs the "smith-waterman" program to do alignment.  Returns 
# an object of type NeedlemanAlignment.
#
# Attributes:
#   executableDir    : i.e., "/home/bmajoros/GF/twain/alignment"
#   matrixFileNoPath : i.e., "blosum62"
# Methods:
#   $aligner=new SmithWaterman($executableDir,$matrixFileNoPath);
#   $alignment=$aligner->align($fasta1,$fasta2,"DNA|PROTEIN",
#                              $gapOpenPenalty,$gapExtendPenalty);
#   $alignment=$aligner->alignSeqs($seq1,$seq2,"DNA|PROTEIN",$open,$extend);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
# $aligner=new SmithWaterman($executableDir,$matrixFileNoPath);
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
#   $alignment=$aligner->alignSeqs($seq1,$seq2,"DNA|PROTEIN",$open,$extend);
sub alignSeqs {
  my ($self,$seq1,$seq2,$type,$gapOpenPenalty,$gapExtendPenalty)=@_;
  my $file1=TempFilename::generate();
  my $file2=TempFilename::generate();
  my $w=new FastaWriter();
  $w->writeFasta(">1",$seq1,$file1);
  $w->writeFasta(">2",$seq2,$file2);
  my $alignment=
    $self->align($file1,$file2,$type,$gapOpenPenalty,$gapExtendPenalty);
  system("rm $file1 $file2");
  return $alignment;
}
#---------------------------------------------------------------------
# $alignment=$aligner->align($fasta1,$fasta2,"DNA",$gapOpen,$gapExtend);
sub align
  {
    my ($self,$fasta1,$fasta2,$type,$gapOpenPenalty,$gapExtendPenalty)=@_;
    my $executableDir=$self->{executableDir};
    my $matrixFileNoPath=$self->{matrixFileNoPath};
    my $exec="$executableDir/smith-waterman";
    my $matrix="$executableDir/matrices/$matrixFileNoPath";
    my $command=
      "$exec $matrix $gapOpenPenalty $gapExtendPenalty $fasta1 $fasta2 $type";
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

