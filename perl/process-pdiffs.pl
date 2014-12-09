#!/usr/bin/perl
#Make SAM-type strings for all pdiffs corresponding to a given gene by calling
#pdiff-to-string.pl

use strict;
use warnings;

die "process-pdiffs.pl <file-prefix> <script-dir>" unless @ARGV==2;

my ($prefix, $scriptDir) = @ARGV;

my $files = `ls $prefix*.pdiff`;
my @files = split(/\s+/, $files);

for my $file (@files)
{
  chomp($file);
  system("perl $scriptDir/pdiff-to-string.pl $file");
}
