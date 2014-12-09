package SamReader;
use strict;
use FileHandle;
use Carp;
use SamRead;

######################################################################
#
# SamReader.pm 
#
# bmajoros@duke.edu 7/11/2011
#
# Reads SAM files.
# 
#
# Attributes:
#   file : FileHandle
#   save : next defline
#   shouldUppercase : bool
# Methods:
#   $reader=new SamReader($filename);
#   $reader=readerFromFileHandle($fileHandle);
#   $read=$reader->nextRead();
#   $reader->close();
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $reader=new SamReader($filename);
sub new
{
  my ($class,$filename)=@_;
  if(!-e $filename) {die $filename." does not exist"}
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
    };
  bless $self,"SamReader";

  return $self;
}
#---------------------------------------------------------------------
#   ($defline,$sequence)=$reader->nextRead();
sub nextRead
  {
    my ($self)=@_;
    my $fh=$self->{file};
    while(!eof($fh)) {
      my $line=<$fh>;
      next if($line=~/^\@/);
      my @fields=split/\t/,$line;
      next unless @fields>=11;
#773_261_1907_F3 16      Chr2    5742023 255     50M     *       0       0       TTATCTGCTCCTACTATTCAATTTCACTCTCAGAGTCGCTGGAAACGCTG      WWP:5DM>:RWTFCOUVFFUSC@OXYPFJE:46KYYTNU^_\GHXZ][ZZ      NM:i:0  NH:i:1
      my ($id,$flags,$substrate,$pos,$mqual,$cigar,$rnext,$pnext,$tlen,$seq,$qual)=@fields;
      --$pos;
      my $read=new SamRead($id,$flags,$substrate,$pos,$mqual,$cigar,
			   $rnext,$pnext,$tlen,$seq,$qual);
      return $read;
    }
    return undef;
  }
#---------------------------------------------------------------------
#   $reader->close();
sub close
  {
    my ($self)=@_;
    close($self->{file});
  }
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

