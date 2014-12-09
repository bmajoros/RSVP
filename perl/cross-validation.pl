#!/usr/bin/perl
#A wrapper script for training-wrapper.pl that performs twofold cross-validation.
#Basically this script calls training-wrapper.pl twice, switching the "chunks" and
#"testchunks" directories, so that the test set becomes the training set (and vice
#versa) for the second run.

use warnings;
use strict;

die "cross-validation.pl <training-param-file>" unless @ARGV==1;

my $paramFile = $ARGV[0];
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

print("Beginning first run: the given chunkdir1 is the training set, and the given chunkdir2 is the test set\n");

#my $firstResult = `perl $scriptDir/training-wrapper.pl $paramFile`;
system("perl $scriptDir/training-wrapper.pl $paramFile");

my $oldChunkDir = $paramHash{"chunkdir1"};
my $oldTestChunkDir = $paramHash{"chunkdir2"};
my $newChunkDir = $oldTestChunkDir;
my $newTestChunkDir = $oldChunkDir;

print("Beginning second run: the given chunkdir2 is the training set, and the given chunkdir1 is the test set\n");

changeParams("chunkdir1", $newChunkDir, "chunkdir2", $newTestChunkDir);

system("perl $scriptDir/training-wrapper.pl $paramFile");

changeParams("chunkdir1", $oldChunkDir, "chunkdir2", $oldTestChunkDir); #reset the original directories

#params: param1, new-value1, param2, new-value2, etc (key-value pairs)
#copied from hillclimb.pl, but using it to change hillclimb parameters here
sub changeParams
{
    my %paramHash = ();
    while(@_ > 0) 
    {
        my $inParam = shift;
        my $inVal = shift;
        $paramHash{$inParam} = $inVal;
    }

    open(OLD_PARAM_FILE, $paramFile);
    my @output = <OLD_PARAM_FILE>;
    close(OLD_PARAM_FILE);
    open(NEW_PARAM_FILE, "> $paramFile");
    for my $line (@output)
    {
        chomp $line;
        my $flag = 0;
        my @keys = keys %paramHash;
        for my $param (@keys)
        {
            my $val = $paramHash{$param};
            if ($flag == 1)
            {
                last;
            }
            if ($line =~ /$param=/)
            {
                print(NEW_PARAM_FILE "$param=$val\n");
                $flag = 1;
            }
        }
        if ($flag == 0)
        {
                print(NEW_PARAM_FILE "$line\n");
        }
    }
}

