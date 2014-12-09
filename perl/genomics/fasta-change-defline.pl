#!/usr/bin/perl
use strict;
use TempFilename;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.fasta> defline ...    (no '>')";
die "$usage\n" unless @ARGV>=2;
#my ($filename,$defline)=@ARGV;
my $filename=shift @ARGV;
my $defline=join(' ',@ARGV);
$defline=">$defline";

my $tempfile=TempFilename::generate();
open(OUT,">$tempfile") || die "Can't create temp file: $tempfile\n";
open(IN,$filename) || die "Can't open $filename\n";
<IN>;
print OUT "$defline\n";
while(<IN>) {print OUT}
close(IN);
close(OUT);

rename($tempfile,$filename) || 
  die "Error renaming $tempfile as $filename\n";


