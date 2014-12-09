package bedfile;
use strict;
use FileHandle;
use bedrecord;

######################################################################
#
# bedfile.pm bmajoros@duke.edu 6/19/2006
#
# 
# 
#
# Attributes:
#
# Methods:
#   $bedfile=new bedfile($filename);
#   my $record=$bedfile->nextRecord(); # returns a "bedrecord" object
#                                      # or undef if end-of-file
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $bedfile=new bedfile($filename);
sub new
{
  my ($class,$filename)=@_;
  
  my $self=
  {
   filehandle=>new FileHandle($filename)
  };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   my $record=$bedfile->nextRecord(); # returns a "bedrecord" object
sub nextRecord
{
    my ($self)=@_;
    my $filehandle=$self->{filehandle};
    my @fields;
    while(@fields!=6)
    {
        my $line=<$filehandle>;
        chomp $line;
        if(!defined($line)) {return undef}
        @fields=split/\t/,$line;
    }
    my ($chrom,$begin,$end,$name,$zero,$strand)=@fields;
    return new bedrecord($chrom,$begin,$end,$name,$strand);
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

