package BlastAligner;
use strict;
use FastaWriter;
use FastaReader;
use TempFilename;

######################################################################
#
# BlastAligner.pm bmajorostigr.org 2/12/2003
#
# Uses BLAST to compute an alignment between two nucleotide sequences.
# 
#
# Attributes:
#
# Methods:
#   $aligner=new BlastAligner();
#   $alignment=$aligner->align(\$seq1,\$seq2,$label1,$label2);
#   $alignment=$aligner->alignFasta($filename1,$filename2,$label1,$label2);
# Private methods:
#   $self->writeFasta(\$seq,$filename,$label);
#   $self->runBlast($file1,$file2,$outfile);
#   my $alignment=$self->parseBlastReport($filename);
#   $self->formatdb($fastaFile);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self=
    {
     fastaWriter=>new FastaWriter,
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $alignment=$aligner->align(\$seq1,\$seq2,$label1,$label2);
sub align
  {
    my ($self,$seq1,$seq2,$label1,$label2)=@_;
    my $file1="$label1.seq";#TempFilename::generate();
    my $file2="$label2.seq";#TempFilename::generate();
    my $outfile="$label1\-$label2.blast";#TempFilename::generate();
    $self->writeFasta($seq1,$file1,$label1);
    $self->writeFasta($seq2,$file2,$label2);
    $self->runBlast($file1,$file2,$outfile);
    my $alignment=$self->parseBlastReport($outfile);
    #unlink $file1;
    #unlink $file2;
    #unlink $outfile;
    return $alignment;
  }
#---------------------------------------------------------------------
#   $alignment=$aligner->alignFasta($filename1,$filename2,$label1,$label2);
sub alignFasta
  {
    my ($self,$filename1,$filename2,$label1,$label2)=@_;

    my $reader=new FastaReader($filename1);
    my ($defline1,$sequence1)=$reader->nextSequence();

    $reader=new FastaReader($filename2);
    my ($defline2,$sequence2)=$reader->nextSequence();

    return $self->align(\$sequence1,\$sequence2,$label1,$label2);
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#   $self->writeFasta(\$seq,$filename,$label);
sub writeFasta
  {
    my ($self,$seq,$filename,$label)=@_;
    my $fastaWriter=$self->{fastaWriter};
    $fastaWriter->writeFastaFromRef(">$label",$seq,$filename);
  }
#---------------------------------------------------------------------
#   $self->formatdb($fastaFile);
sub formatdb
  {
    my ($self,$filename)=@_;
    my $command="formatdb -i $filename -p F";
    return system($command);
  }
#---------------------------------------------------------------------
#   $self->runBlast($file1,$file2,$outfile);
sub runBlast
  {
    my ($self,$file1,$file2,$outfile)=@_;
    $self->formatdb($file1);
    my $command="blastall -p blastn -d $file1 -i $file2 > $outfile";
    return system($command);
  }
#---------------------------------------------------------------------
#   my $alignment=$self->parseBlastReport($filename);
sub parseBlastReport
  {
    my ($self,$filename)=@_;
    open(IN,$filename) || die "Can't open file: \"$filename\"";
    my $ignore=1;
    my @alignment;
    while(<IN>)
      {
	if(/Strand\s*=\s*(\S+)\s*\/\s*(\S+)/)
	  {
	    my ($s1,$s2)=("\L$1","\L$2");
	    $ignore=($s1 ne "plus" || $s2 ne "plus");
	  }
	elsif(/^Query:\s*(\d+)\s*(\S+)\s*(\d+)\s*$/)
	  {
	    next if $ignore;
	    my ($queryBegin,$querySeq,$end)=($1-1,$2,$3);
	
	    $_=<IN>;
	    $_=~/^.{12}(.*)/;
	    my $alignString=$1;

	    $_=<IN>;
	    $_=~/^Sbjct:\s*(\d+)\s*(\S*)\s*\d+\s*$/;
	    my ($subjectBegin,$subjectSeq)=($1-1,$2);
	
	    #print "$querySeq\n$alignString\n$subjectSeq\n";

	    my $queryPos=$queryBegin;
	    my $subjectPos=$subjectBegin;
	    my $subjectStringLength=length $subjectSeq;
	    for(my $i=0 ; $i<$subjectStringLength ; ++$i)
	      {
		my $s=substr($subjectSeq,$i,1);
		my $q=substr($querySeq,$i,1);
		if($s ne "-")
		  {
		    $alignment[$subjectPos]=$queryPos;
		    ++$subjectPos;
		  }
		if($q ne "-") {++$queryPos}
	      }
	  }
      }
    close(IN);
    return \@alignment;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

1;

