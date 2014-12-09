package ConfigFile;
use strict;

######################################################################
#
# ConfigFile.pm bmajorostigr.org 5/17/2004
#
# 
# 
#
# Attributes:
#
# Methods:
#   $configFile=new ConfigFile($filename);
#   $value=$configFile->lookup($key);
#   $value=$configFile->lookupOrDie($key);
#   
# Private Methods:
#   $self->load($filename);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $configFile=new ConfigFile($filename);
sub new
{
  my ($class,$filename)=@_;
  
  my $self=
    {
     hash=>{},
    };
  bless $self,$class;
  
  $self->load($filename);

  return $self;
}
#---------------------------------------------------------------------
#   $value=$configFile->lookup($key);
sub lookup
  {
    my ($self,$key)=@_;
    return $self->{hash}->{$key};
  }
#---------------------------------------------------------------------
#   $value=$configFile->lookupOrDie($key);
sub lookupOrDie
  {
    my ($self,$key)=@_;
    die "$key not defined in config file\n" unless defined $self->{hash}->{$key};
    return $self->{hash}->{$key};
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $self->load($filename);
sub load
  {
    my ($self,$filename)=@_;
    my $hash=$self->{hash};
    open(IN,$filename) || die "can't open $filename\n";
    while(<IN>)
      {
	$_=~s/\#.*//g;
	if(/\s*(\S+)\s*=\s*(\S+)/)
	  {
	    $hash->{$1}=$2;
	  }
      }
    close(IN);
  }


1;

