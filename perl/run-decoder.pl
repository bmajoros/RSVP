#!/usr/bin/perl
use warnings;
use strict;

die "run-decoder.pl <data-prefix> <out-prefix> <param-file>" unless @ARGV == 3;

my $dataPrefix = $ARGV[0];
my $outPrefix = $ARGV[1];
my $paramFile = $ARGV[2];

open(PARAM, $paramFile);
my %paramHash = ();

while(<PARAM>)
{
    chomp;
    if (/^\s*$/) {next;}
    if (/^#/) {next;}
    my ($key, $val) = split(/=/);
    if (!defined($val) || !defined($key)) {next;}
    $paramHash{$key} = $val;
}

my $scriptDir = $paramHash{"scriptdir"};
my $javaDir = $paramHash{"classdir"};

system("gunzip $dataPrefix.graph.gz");
system("java -Xmx512m -classpath $javaDir RunDecoder $dataPrefix $paramFile > $outPrefix.gff");
system("gzip $dataPrefix.graph");
