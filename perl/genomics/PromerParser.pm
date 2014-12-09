package PromerParser;
use strict;
use PromerMatch;

######################################################################
#
# PromerParser.pm bmajorostigr.org 8/22/2003
#
# 
# 
#
# Attributes:
#
# Methods:
#   $pp=new PromerParser();
#   $matches=$pp->parse($filename);
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
#   $matches=$pp->parse($filename);
sub parse
  {
    my ($self,$filename)=@_;
    open(IN,$filename) || 
      die "Can't open $filename in PromerParser::parse()";
    my $seenEquals=0;
    my $matches=[];
    while(<IN>)
      {
	if(/===/) {$seenEquals=1}
	elsif($seenEquals)
	  {
	    $_=~/(\S.*\S)/ || next;
	    $_=$1;
	    my @fields=split/\s+/,$_;
	    my ($begin1,$end1,$d1,$begin2,$end2,$d2,$len1,$len2,
		$d3,$pctIdentity,$pctSimilarity,$pctStop,$d4,
		$frame1,$frame2,$id1,$id2)=@fields;
	    my $strand1='+';
	    my $strand2='+';
	    if($begin1>$end1)
	      {
		#($begin1,$end1)=($end1,$begin1);
		$strand1='-';
	      }
	    if($begin2>$end2)
	      {
		#($begin2,$end2)=($end2,$begin2);
		$strand2='-';
	      }

	    my $match=new PromerMatch($begin1,$end1,$begin2,$end2,
				      $len1,$len2,$strand1,$strand2,
				      $pctIdentity,$pctSimilarity,
				      $pctStop,$id1,$id2,$frame1,
				      $frame2);
	    push @$matches,$match;
	  }
      }
    close(IN);
    return $matches;
  }
#---------------------------------------------------------------------
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

1;

