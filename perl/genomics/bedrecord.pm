package bedrecord;
use strict;

######################################################################
#
# bedrecord.pm bmajoros@duke.edu 6/19/2006
#
#   Represents a single record from a BED file.  See bedfile.pm.
# 
#
# Attributes:
#   chrom
#   begin
#   end
#   name
#   strand
# Methods:
#   $r=new bedrecord($chrom,$begin,$end,$name,$strand);
#   $chrom=$r->getChrom();
#   $begin=$r->getBegin();
#   $end=$r->getEnd();
#   $name=$r->getName();
#   $strand=$r->getStrand();
#   $line=$r->asString();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $r=new bedrecord($chrom,$begin,$end,$name,$strand);
sub new
{
  my ($class,$chrom,$begin,$end,$name,$strand)=@_;

  my $self=
  {
   chrom=>$chrom,
   begin=>$begin,
   end=>$end,
   name=>$name,
   strand=>$strand
  };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $chrom=$r->getChrom();
sub getChrom
{
    my ($self)=@_;
    return $self->{chrom};
}
#---------------------------------------------------------------------
#   $begin=$r->getBegin();
sub getBegin
{
    my ($self)=@_;
    return $self->{begin};
}
#---------------------------------------------------------------------
#   $end=$r->getEnd();
sub getEnd
{
    my ($self)=@_;
    return $self->{end};
}
#---------------------------------------------------------------------
#   $name=$r->getName();
sub getName
{
    my ($self)=@_;
    return $self->{name};
}
#---------------------------------------------------------------------
#   $strand=$r->getStrand();
sub getStrand
{
    my ($self)=@_;
    return $self->{strand};
}
#---------------------------------------------------------------------
#   $line=$r->asString();
sub asString
{
    my ($self)=@_;
    return "$self->{chrom}\t$self->{begin}\t$self->{end}\t$self->{name}\t0\t$self->{strand}";
}
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

