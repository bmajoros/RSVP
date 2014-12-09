#!/usr/bin/perl
use strict;
use ProgramName;

# bmajoros 12/11/12

my $name=ProgramName::get();
die "$name <chunks-dir> <chrom-pileup-dir> <package-dir>\n" unless @ARGV==3;
my ($chunksDir,$pileupDir,$packageDir)=@ARGV;

my @files=`ls $chunksDir`;
my $n=@files;
for(my $i=0 ; $i<$n; ++$i) {
  my $file=$files[$i];
  chomp $file;
  next unless $file=~/([^\/]+)\.coords\.gff/;
  my $chunkID=$1;
  my ($substrate,$begin,$end,$strand);
  open(IN,"$chunksDir/$file") || die;
  while(<IN>) {
    chomp;
    my @fields=split/\t/,$_;
    $substrate=$fields[0];
    $begin=$fields[3]-1;
    $end=$fields[4];
    $strand=$fields[6];
  }
  close(IN);
  my $cmd="$packageDir/c++/extract-pileup-interval $pileupDir/$substrate.pileup $pileupDir/$substrate.junctions $begin $end $chunksDir/$chunkID";
  #print "$cmd\n\n";
  system($cmd);
}

