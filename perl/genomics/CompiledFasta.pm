package CompiledFasta;
use strict;
use FileHandle;

######################################################################
#
# CompiledFasta.pm bmajorostigr.org 2/24/2005
#
# For reading a FASTA file that has no defline or whitespace; only
# sequence.  This allows us to extract substrings by indexing into
# the file using fseek().
#
# Attributes:
#   file
# Methods:
#   $compiledFasta=new CompiledFasta($filename);
#   $str=$compiledFasta->load($begin,$len);
#   $compiledFasta->close();
#   $size=$compiledFasta->fileSize();
#   
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $compiledFasta=new CompiledFasta($filename);
sub new
{
  my ($class,$filename)=@_;
  
  my $self=
    {
     file=>new FileHandle($filename),
     name=>$filename
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $str=$compiledFasta->load($begin,$len);
sub load
  {
    my ($self,$begin,$len)=@_;
    my $file=$self->{file};
    seek($file,$begin,0);
    my $buffer;
    read($file,$buffer,$len);
    return $buffer;
  }
#---------------------------------------------------------------------
#   $compiledFasta->close();
sub close
  {
    my ($self)=@_;
    $self->{file}=0;
    undef $self->{file};
  }
#---------------------------------------------------------------------
sub fileSize
  {
    my ($self)=@_;
    my $filename=$self->{name};
    my @stats=stat $filename;
    return $stats[7];
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

