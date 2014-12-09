#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <input-chunks> <out-multi> <out-single>\n" unless @ARGV==3;
my ($inputChunks,$multiChunks,$singleChunks)=@ARGV;

system("mkdir -p $multiChunks $singleChunks");

my $reader=new GffTranscriptReader;
my @files=`ls $inputChunks`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless $file=~/^(\d+)\.gff$/;
  my $id=$1;
  my $path="$inputChunks/$file";
  my $transcripts=$reader->loadGFF($path);
  #my $geneList=$reader->loadGenes($path);
  #die unless $geneList;
  #print "$path\n";
  #die @$geneList unless @$geneList==1;
  #my $gene=$geneList->[0];
  #my $numTranscripts=$gene->getNumTranscripts();
  my $numTranscripts=@$transcripts;
  die unless $numTranscripts>0;
  my $outDir=$numTranscripts==1 ? $singleChunks : $multiChunks;
  system("cp $inputChunks/$id.* $outDir");
}


