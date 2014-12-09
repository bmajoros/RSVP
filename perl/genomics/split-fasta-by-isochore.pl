
unfinished


#!/usr/bin/perl
use strict;
use ProgramName;
use GffReader;
use FastaWriter;

# Process command line
my $name=ProgramName::get();
my $usage="
$name <*.fasta> <isochores.gff> <max-chunk-size>

  NOTE: the fasta file must contain only one sequence, and the
        GFF features must all occur on that substrate
";
die "$usage\n" unless @ARGV==3;
my ($fastaFile,$isochoreFile,$maxChunkSize)=@ARGV;
if($fastaFile=~/\.cof$/) 
  {die "compiled FASTA files are not yet supported\n"}

# Load isochore boundaries
my $gffReader=new GffReader;
my $isochores=$gffReader->loadGFF();
my $numIsochores=@$isochores;

# Process input FASTA file
my $writer=new FastaWriter;
my $currentIsoIndex=0;
my $currentIsochore=$isochores->[$currentIsoIndex];
my $nextBoundary=$currentIsochore->getEnd();
my $offset=0;
my $buffer="";
my $chunkId=1;
open(IN,$fastaFile) || die "Can't open $fastaFile\n";
my $defline=<IN>;
while(<IN>)
  {
    


  }
close(IN);



