#!/usr/bin/perl
use strict;
use ProgramName;
use Getopt::Std;

our ($opt_o,$opt_c);
getopts("oc");

my $name=ProgramName::get();
die "
$name [-oc] <infile> <filestem-prefix>

where:
       -o = overwrite files even if they exist
       -c = compress output files
" 
unless @ARGV==2;
my ($infile,$prefix)=@ARGV;

my $isOpen=0;
open(IN,$infile) || die "can't open file $infile\n";
my @filenames;
while(<IN>)
{
    if(/^>(\S+)/)
    {
        my $id=$1;
        if($isOpen) {close(OUT)}
        my $filename="$prefix$id.fasta";
        if(-e $filename && !$opt_o) 
          {die "$filename exists; use -o to overwrite\n"}
        open(OUT,">$filename") || die "can't create file $filename\n";
        push @filenames,$filename;
        $isOpen=1;
        print OUT $_;
    }
    else 
    {
        print OUT $_;
    }
}
close(IN);
if($isOpen) {close(OUT)}

if($opt_c)
{
    foreach my $infile (@filenames)
    {
        $infile=~/(.*)\.fasta/ || die;
        my $outfile="$1.dna";
        system("compress-fasta $infile $outfile");
        system("rm $infile");
    }
}



