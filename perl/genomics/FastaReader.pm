package FastaReader;
use strict;
use FileHandle;
use Carp;

######################################################################
#
# FastaReader.pm 
#
# bmajoros@tigr.org 7/8/2002
#
# Reads fasta files.
#
#
# Attributes:
#   file : FileHandle
#   save : next defline
#   shouldUppercase : bool
# Methods:
#   $reader=new FastaReader($filename);
#   $reader=readerFromFileHandle($fileHandle);
#   ($defline,$sequence)=$reader->nextSequence();
#   ($defline,$seqRef)=$reader->nextSequenceRef();
#   $reader->close();
#   $size=FastaReader::getGenomicSize($filename);
#   FastaReader::readAll($filename); # returns hash : id->sequence
#   FastaReader::readAllAndKeepDefs($filename); # returns hash : id->[def,seq]
#   $reader->dontUppercase();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $reader=new FastaReader($filename);
sub new
{
  my ($class,$filename)=@_;
  if(!-e $filename) {die "$filename does not exist"}
  return readerFromFileHandle(new FileHandle($filename));
}
#---------------------------------------------------------------------
#   $reader=readerFromFileHandle($fileHandle);
sub readerFromFileHandle
{
  my ($fileHandle)=@_;
  if(!$fileHandle) {die "readerFromFileHandle(NULL)"}
  
  my $self=
    {
     file=>$fileHandle,
     shouldUppercase=>1,
    };
  bless $self,"FastaReader";

  return $self;
}
#---------------------------------------------------------------------
#   ($defline,$sequence)=$reader->nextSequence();
sub nextSequence
  {
    my ($self)=@_;
    my ($defline,$seqRef)=$self->nextSequenceRef();
    return $seqRef ? ($defline,$$seqRef) : ($defline,$seqRef);
  }
#---------------------------------------------------------------------
#   ($defline,$seqRef)=$reader->nextSequenceRef();
sub nextSequenceRef
  {
    my ($self)=@_;

    die unless defined $self;
    my $file=$self->{file};
    return (undef,undef) if(eof($file));

    my $defline="";
    my $seq="";
    while(defined($self->{save}) || ($_=<$file>))
      {
	$_=$self->{save} if defined($self->{save});
	if(/^\>/)
	  {
	    $defline=$_;
	    while(<$file>)
	      {
		if(/^\>/)
		  {
		    $self->{save}=$_;
		    last;
		  }
		$_=~s/\s+//g;
                if($self->{shouldUppercase}) {$seq.="\U$_"}
                else {$seq.="$_"}
	      }
	    return ($defline,\$seq);
	  }
      }
    return (undef,undef);
  }
#---------------------------------------------------------------------
#   $reader->close();
sub close
  {
    my ($self)=@_;
    close($self->{file});
  }
#---------------------------------------------------------------------
#   $size=FastaReader::getGenomicSize($filename);
sub getGenomicSize
  {
    my ($filename)=@_;

    my $size=0;
    my $reader=new FastaReader($filename);
    while(1)
      {
	my ($defline,$sequence)=$reader->nextSequence();
	last unless $defline;
	$size+=length($sequence);
      }
    return $size;
  }
#---------------------------------------------------------------------
sub readAll
  {
    my ($filename)=@_;
    my $reader=new FastaReader($filename);
    my %hash;
    while(1)
      {
	my ($def,$seq)=$reader->nextSequence();
	last unless defined $def;
	$def=~/>\s*(\S+)/;
	my $id=$1;
	if($seq=~/\n/) {chop $seq}
	$hash{$id}=$seq;
      }
    return \%hash;
  }
#---------------------------------------------------------------------
#   FastaReader::readAllAndKeepDefs($filename); # returns hash : id->[def,seq]
sub readAllAndKeepDefs
  {
    my ($filename)=@_;
    my $reader=new FastaReader($filename);
    my %hash;
    while(1)
      {
	my ($def,$seq)=$reader->nextSequence();
	last unless defined $def;
	$def=~/>\s*(\S+)/;
	my $id=$1;
	if($seq=~/\n/) {chop $seq}
	$hash{$id}=[$def,$seq];
      }
    return \%hash;
  }
#---------------------------------------------------------------------
#   $reader->dontUppercase();
sub dontUppercase
  {
    my ($self)=@_;
    $self->{shouldUppercase}=0;
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

