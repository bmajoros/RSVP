#!/usr/bin/perl
#This script runs hillclimb.pl repeatedly, with random initial values
#for the parameters. The number of times to repeat the hill-climbing
#program, as well as the range of possible lambda and gamma values,
#are stored in the training parameter file (the same file used by
#hillclimb.pl). Uses chunkdir1 as the training set and chunkdir2
#as the test set.
#
#Currently, lambda3's value is hardcoded on line 52, since it is not
#trained by hillclimb.pl. Edit this line to change its value.

use warnings;
use strict;

die "training-wrapper.pl <training-param-file>" unless @ARGV == 1;

my $prog = "hillclimb.pl";
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

close(PARAM);

my $scriptDir = $paramHash{"scriptdir"};
my $lambdaRange = $paramHash{"lambdarange"};
my $gammaRange = $paramHash{"gammarange"};
my $gammaMin = $paramHash{"gammamin"};
my $numRuns = $paramHash{"numruns"};

my ($bestLambda2, $bestLambda3, $bestBeta, $bestGamma) = (0,0,0,0);
my $bestFscore = 0;

for (my $i = 0; $i < $numRuns; $i++)
{
    my $startTime = `date`;
    chomp($startTime);
    my $runNum = $i + 1;
    print("run $runNum started at $startTime\n");

    my $lambda2 = rand($lambdaRange);
    my $lambda3 = 1;
    my $beta = rand(1);
    my $gamma = rand($gammaRange) + $gammaMin;
    print("starting values: lambda2 = $lambda2, lambda3 = $lambda3, beta = $beta, gamma = $gamma\n");
    changeParams("startlambda2", $lambda2, "startlambda3", $lambda3, "startbeta", $beta, "startgamma", $gamma);
    
    my $result = `perl $scriptDir/$prog $paramFile 2>/dev/null  | tail -n 1`;
    chomp($result);
    print("$result\n");

    my $endTime = `date`;
    chomp($endTime);
    print("run $runNum ended at $endTime\n");

    my %resultHash = ();
    my @resultSplit = split(/\s+/, $result);
    for my $result (@resultSplit)
    {
        my ($key, $val) = split(/=/, $result);
        $resultHash{$key} = $val;
    }
    my $fscore = $resultHash{"f"};
    if ($fscore > $bestFscore)
    {
        $bestLambda2 = $resultHash{"lambda2"};
        $bestLambda3 = $resultHash{"lambda3"};
        $bestBeta = $resultHash{"beta"};
        $bestGamma = $resultHash{"gamma"};
        $bestFscore = $fscore;
    }
}

my $chunkDir = $paramHash{"chunkdir1"};
my $testChunkDir = $paramHash{"chunkdir2"};
my $N = $paramHash{"N"};
changeParams("startlambda2", $bestLambda2, "startlambda3", $bestLambda3, "startbeta", $bestBeta,
    "startgamma", $bestGamma, "chunkdir1", $testChunkDir, "N", 0);
my $testResult = `perl $scriptDir/$prog $paramFile 2>/dev/null  | tail -n 1`;
chomp($testResult);
print("test result: $testResult\n");

changeParams("N", $N, "chunkdir1", $chunkDir); #reset params to old values--may be unecessary

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
            if ($line =~ /^$param=/)
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

