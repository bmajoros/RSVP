#!/usr/bin/perl
use strict;
use ProgramName;
use FastaWriter;
use Getopt::Std;
use CompiledFasta;
our $opt_c;
getopts('c:');

my $name=ProgramName::get();
my $usage="
$name [-c ID] <fasta> <chunk-size> <out-dir> <mapping-file.out>

  where -c indicates that <fasta> is compiled (*.cof), substrate=ID
           (compiled fasta files contain no defline and no whitespace)
";
die "$usage\n" unless @ARGV==4;
my ($infile,$chunkSize,$outDir,$mappingFile)=@ARGV;

my $fastaWriter=new FastaWriter;
open(MAP,">$mappingFile") || die "can't create file: $mappingFile\n";
print MAP "chunk\tcontig\toffset\n";
my $isCompiled=$opt_c;
if(!$isCompiled && $infile=~/\.cof$/) 
  {die "\nuse -c option with compiled fastas!\n\n"}

#################### COMPILED FASTA FILES ########################
if($isCompiled)
  {
    my $fasta=new CompiledFasta($infile);
    my $fileSize=$fasta->fileSize();
    my $substrateId=$opt_c;
    my $chunkId=1;
    my $offset=0;
    while($offset+$chunkSize<=$fileSize)
      {
	my $chunk=$fasta->load($offset,$chunkSize);
	my $outfile="$outDir/${substrateId}_$chunkId.cof";
	open(OUT,">$outfile") || die "can't create file: $outfile\n";
	print OUT $chunk;
	close(OUT);
	print MAP "${substrateId}_$chunkId\t$substrateId\t$offset\n";
	$offset+=$chunkSize;
	++$chunkId;
      }
    my $finalChunkSize=$fileSize-$offset;
    if($finalChunkSize>0)
      {
	my $chunk=$fasta->load($offset,$finalChunkSize);
	my $outfile="$outDir/${substrateId}_$chunkId.cof";
	open(OUT,">$outfile") || die "can't create file: $outfile\n";
	print OUT $chunk;
	close(OUT);
	print MAP "${substrateId}_$chunkId\t$substrateId\t$offset\n";
      }
  }

################## NON-COMPILED FASTA FILES #####################
else
  {
    open(IN,$infile) || die "can't open file: $infile\n";
    my $defline=<IN>;
    $defline=~/^>\s*(\S+)/ || die "can't parse defline in file $infile\n";
    my $substrateId=$1;
    my $chunkId=1;
    my $offset=0;
    my ($chunk,$len,$buffer);
    while(<IN>)
      {
	$_=~s/\s+//g;
	$buffer.=$_;
	my $maxRoom=$chunkSize-$len;
	my $bufferLen=length $buffer;
	if($bufferLen<$maxRoom) 
	  {$chunk.=$buffer;$len+=$bufferLen;$buffer="";next}
	my $piece=substr($buffer,0,$maxRoom);
	$chunk.=$piece;
	$buffer=substr($buffer,$maxRoom,$bufferLen-$maxRoom);
	my $defline=
	  ">${substrateId}_$chunkId /contig=$substrateId /offset=$offset";
	my $outfile="$outDir/${substrateId}_$chunkId.fasta";
	$fastaWriter->writeFasta($defline,$chunk,$outfile);
	print MAP "${substrateId}_$chunkId\t$substrateId\t$offset\n";
	$offset+=length $chunk;
	$len=0;
	$chunk="";
	++$chunkId;
      }
    my $defline=
      ">${substrateId}_$chunkId /contig=$substrateId /offset=$offset";
    my $outfile="$outDir/${substrateId}_$chunkId.fasta";
    $fastaWriter->writeFasta($defline,$chunk,$outfile);
    print MAP "${substrateId}_$chunkId\t$substrateId\t$offset\n";
    close(IN);
  }
close(MAP);
