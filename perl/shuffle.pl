#!/usr/bin/perl
#shuffle the lines of a given input file

use warnings;
use strict;
use List::Util 'shuffle';

my @input;

if (@ARGV==1)
{
    my $inputFileName = shift;
    open(IN_FILE, $inputFileName);
    @input = <IN_FILE>;
}
else
{
    @input = <STDIN>;
}

my @output = shuffle(@input);

print("@output");
