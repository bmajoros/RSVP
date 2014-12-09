#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in-dir> <out-dir> <n>\n" unless @ARGV==3;
my ($inDir,$outDir,$n)=@ARGV;

my @files=`ls $inDir`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  process($files[$i]);
}

sub process {
  my ($filename)=@_;
  chomp $filename;
  my $infile="$inDir/$filename";
  my $outfile="$outDir/$filename";
  open(IN,$infile) || die "Can't open file: $infile\n";
  open(OUT,">$outfile") || die "Can't create file: $outfile\n";
  while(<IN>) {
    if(/parse=(\d+)/) {
      my $index=$1;
      if($index<=$n) { print OUT }
    }
    elsif(/transcript_id=\"CUFF.\d+.(\d+)\";/) {
      my $index=$1;
      my @fields=split/\t/,$_;
      if($fields[3]<0) {$fields[3]=0}
      $_=join("\t",@fields);
      if($index<=$n) { print OUT }
    }
  }
  close(OUT);
  close(IN);
}




