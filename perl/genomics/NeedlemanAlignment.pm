package NeedlemanAlignment;
use strict;

######################################################################
#
# NeedlemanAlignment.pm bmajoros@tigr.org 5/17/2004
#
# 
# 
#
# Attributes:
#   alignment=>array[3] where [0] & [1] are gapped sequences and
#                             [2] contains spaces and '|' characters
#   file1=>filename of first sequence
#   file2=>filename of second sequence
#   percentIdentity=>integer
#   matches=>integer
#   mismatches=>integer
#   insertions=>integer
#   length=>integer
#   score=>integer
#
# Methods:
#   $aln=new NeedlemanAlignment($alignment,$seq1name,$seq2name,
#                               $percentIdentity,$matches,$mismatches,
#                               $insertions,$length,$score);
#   $alignment=$aln->getAlignment();
#   $seq1name=$aln->getSeq1Name();
#   $seq2name=$aln->getSeq2Name();
#   $pctIdent=$aln->getPercentIdentity();
#   $matches=$aln->getNumMatches();
#   $mismatches=$aln->getNumMismatches();
#   $insertions=$aln->getNumInsertions();
#   $len=$aln->getLength();
#   $score=$aln->getScore();
#   $aln->print($pageWidth);
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $aln=new NeedlemanAlignment($alignment,$seq1name,$seq2name,
#                               $percentIdentity,$matches,$mismatches,
#                               $insertions,$length,$score);
sub new
{
  my ($class,$alignment,$seq1name,$seq2name,
      $percentIdentity,$matches,$mismatches,
      $insertions,$length,$score)=@_;
  
  my $self=
    {
     alignment=>$alignment,
     file1=>$seq1name,
     file2=>$seq2name,
     percentIdentity=>$percentIdentity,
     matches=>$matches,
     mismatches=>$mismatches,
     insertions=>$insertions,
     length=>$length,
     score=>$score
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $alignment=$aln->getAlignment();
sub getAlignment
  {
    my ($self)=@_;
    return $self->{alignment};
  }
#---------------------------------------------------------------------
#   $seq1name=$aln->getSeq1Name();
sub getSeq1Name
  {
    my ($self)=@_;
    return $self->{file1};
  }
#---------------------------------------------------------------------
#   $seq2name=$aln->getSeq2Name();
sub getSeq2Name
  {
    my ($self)=@_;
    return $self->{file2};
  }
#---------------------------------------------------------------------
#   $pctIdent=$aln->getPercentIdentity();
sub getPercentIdentity
  {
    my ($self)=@_;
    return $self->{percentIdentity};
  }
#---------------------------------------------------------------------
#   $matches=$aln->getNumMatches();
sub getNumMatches
  {
    my ($self)=@_;
    return $self->{matches};
  }
#---------------------------------------------------------------------
#   $mismatches=$aln->getNumMismatches();
sub getNumMismatches
  {
    my ($self)=@_;
    return $self->{mismatches};
  }
#---------------------------------------------------------------------
#   $insertions=$aln->getNumInsertions();
sub getNumInsertions
  {
    my ($self)=@_;
    return $self->{insertions};
  }
#---------------------------------------------------------------------
#   $len=$aln->getLength();
sub getLength
  {
    my ($self)=@_;
    return $self->{length};
  }
#---------------------------------------------------------------------
#   $score=$aln->getScore();
sub getScore
  {
    my ($self)=@_;
    return $self->{score};
  }
#---------------------------------------------------------------------
#   $aln->print($pageWidth);
sub print
  {
    my ($self,$pageWidth)=@_;
    my $alignment=$self->{alignment};
    my $len=length($alignment->[0]);
    my $begin=0;
    while($begin<$len)
      {
	my $top=substr($alignment->[0],$begin,$pageWidth);
	my $middle=substr($alignment->[1],$begin,$pageWidth);
	my $bottom=substr($alignment->[2],$begin,$pageWidth);
	print "Query: $top\n       $middle\nSbjct: $bottom\n\n";
	$begin+=$pageWidth;
      }
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

