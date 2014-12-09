package FastaWriter;
use strict;
use Carp;
use FileHandle;

######################################################################
#
# FastaWriter.pm
#
# bmajoros@tigr.org 7/8/2002
#
# Given a defline and a sequence, formats the sequence to 60-char
# lines and writes into fasta file.
#
# Attributes:
#
# Methods:
#   $w=new FastaWriter($optionalWidth);
#   $w->writeFasta($defline,$sequence,$filename);
#   $w->writeFastaFromRef($defline,\$sequence,$filename);
#   $w->addToFasta($defline,$sequence,$filehandle);
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $w=new FastaWriter($optionalWidth);
sub new
{
  my ($class,$optionalWidth)=@_;
  
  $optionalWidth=60 unless defined $optionalWidth;
  my $self=
    {
     width => $optionalWidth,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $w->writeFasta($defline,$sequence,$filename);
sub writeFasta
  {
      my ($self,$defline,$sequence,$filename)=@_;
      my $h=new FileHandle(">$filename") || die "Can't create $filename";
      $self->addToFasta($defline,$sequence,$h);
      close($h);
  }
#---------------------------------------------------------------------
#   $w->writeFastaFromRef($defline,\$sequence,$filename);
sub writeFastaFromRef
  {
      my ($self,$defline,$seqRef,$filename)=@_;
      my $h=new FileHandle(">$filename") || die "Can't create $filename";
      $self->addToFasta($defline,$$seqRef,$h);
      close($h);
  }
#---------------------------------------------------------------------
#   $w->addToFasta($defline,$sequence,$filehandle);
sub addToFasta
{
    my ($self,$defline,$sequence,$filehandle)=@_;
    
    chop $defline if($defline=~/\n$/);
    $defline=">$defline" unless $defline=~/^>/;

    my $width=$self->{width};
    print $filehandle "$defline\n";
    my $length=length $sequence;
    my $numLines=int($length/$width);
    ++$numLines if($length % $width);
    my $start=0;
    for(my $i=0 ; $i<$numLines ; ++$i)
      {
	my $line=substr($sequence,$start,$width);
	print $filehandle "$line\n";
	$start+=$width;
      }
    if($length==0) {print $filehandle "\n"}
  }
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

