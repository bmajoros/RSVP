#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <chunks1> <chunks2> <outdir>\n" unless @ARGV==3;
my ($dir1,$dir2,$outdir)=@ARGV;

my @files=`ls $dir1`;
my $N=@files;
my $largest=0;
for(my $i=0 ; $i<$N ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  system("cp $dir1/$file $outdir");
  $file=~/(\d+)\.fasta/ || next;
  if($1>$largest) { $largest=$1 }
}

my @files=`ls $dir2`;
my $N=@files;
for(my $i=0 ; $i<$N ; ++$i) {
  my $file=$files[$i];
  chomp $file;
  $file=~/(\d+)\.(\S+)/ || die;
  my $oldID=$1;
  my $newID=$oldID+$largest;
  my $ext=$2;
  system("cp $dir2/$file $outdir/$newID.$ext");
  if($ext eq "coords.gff")
    { system("sub $outdir/$newID.$ext chunk_id=$oldID chunk_id=$newID") }
  elsif($ext eq "fasta")
    { system("sub $outdir/$newID.$ext '>$oldID' '>$newID'") }
}
print "delta=$largest\n";




