#!/usr/bin/perl
#This script is a basic training program, using hill-climbing, for
#the Decoder pipeline. Right now the script just trains the parameter vector
#by optimizing each individual parameter, using its current step size,
#in turn during each iteration. A more advanced training algorithm
#(e.g. RPROP) may be worth using. The changeParams and getFScore subroutines
#should be flexible enough to allow any changes to the training algorithm.
#All training parameters (start values, step sizes, number of iterations,
#and directories) are stored in a training parameter file.
#
#NOTES:
#gamma upper bound is hard-coded in right now (along with other bounds)
#uses chunkdir1 for training
#the script now ignores lambda3 during training

use strict;
use warnings;

die ("hillclimb.pl <training-param-file>") unless @ARGV == 1;

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

my $chunkDir = $paramHash{"chunkdir1"};
my $outBaseDir = $paramHash{"outbasedir"};
my $scriptDir = $paramHash{"scriptdir"};
my $decoderParamFile = $paramHash{"decoderparamfile"};

my $dl2 = $paramHash{"steplambda2"};
my $dl3 = $paramHash{"steplambda3"};
my $db = $paramHash{"stepbeta"};
my $dg = $paramHash{"stepgamma"};

my $lambda2 = $paramHash{"startlambda2"};
my $lambda3 = $paramHash{"startlambda3"};
my $beta = $paramHash{"startbeta"};
my $gamma = $paramHash{"startgamma"};

my $N = $paramHash{"N"};

my %dpHash = ();

my $currentRun = 1;

for (my $i = 1; $i <= $N; $i++) #based on the simple algorithm in Bill's book (pg 41)
{
    my $sl2;
    if (f($lambda2 + $dl2, $lambda3, $beta, $gamma) > f($lambda2, $lambda3, $beta, $gamma))
    {
        $sl2 = $dl2;
    }
    else
    {
        $sl2 = -$dl2;
    }
    while (f($lambda2 + $sl2, $lambda3, $beta, $gamma) >= f($lambda2, $lambda3, $beta, $gamma))
    {
        $lambda2 = $lambda2 + $sl2;
    }
    $dl2 = $dl2 / 2;

    my $sb;
    if (f($lambda2, $lambda3, $beta + $db, $gamma) > f($lambda2, $lambda3, $beta, $gamma))
    {
        $sb = $db;
    }
    else
    {
        $sb = -$db;
    }
    while (f($lambda2, $lambda3, $beta + $sb, $gamma) >= f($lambda2, $lambda3, $beta, $gamma))
    {
        $beta = $beta + $sb;
    }
    $db = $db / 2;

    my $sg;
    if (f($lambda2, $lambda3, $beta, $gamma + $dg) > f($lambda2, $lambda3, $beta, $gamma))
    {
        $sg = $dg;
    }
    else
    {
        $sg = -$dg;
    }
    while (f($lambda2, $lambda3, $beta, $gamma + $sg) >= f($lambda2, $lambda3, $beta, $gamma))
    {
        $gamma = $gamma + $sg;
    }
    $dg = $dg / 2;


}

my $finalf = f($lambda2, $lambda3, $beta, $gamma);

print("f=$finalf lambda2=$lambda2 lambda3=$lambda3 beta=$beta gamma=$gamma\n");

#objective function
#params: lambda2, lambda3, beta, gamma
sub f
{
    my $lambda2 = shift;
    my $lambda3 = shift;
    my $beta = shift;
    my $gamma = shift;


    my $key = "$lambda2 $lambda3 $beta $gamma";

    if ($lambda2 < 0 || $lambda3 < 0 || ($beta < 0 || $beta > 1) || ($gamma < 1 || $gamma > 10))
    {
        print("one or more params out of bounds: $key\n");
        return 0;
    }

    if (exists($dpHash{$key}))
    {
        my $fscore = $dpHash{$key};
        print("used hash for params $key to get fscore $fscore\n");
        return $fscore;
    }

    system("mkdir $outBaseDir/$currentRun");    
    
    changeParams("lambda2", $lambda2, "lambda3", $lambda3, "beta", $beta, "gamma", $gamma, "outdir", "$outBaseDir/$currentRun", "chunkdir", $chunkDir);
    system("perl $scriptDir/batch-decoder-revcomp.pl $decoderParamFile");
    print("current params $key\n");
    print("waiting for SGE jobs to finish...\n");
    while (`qstat` =~ /\w/)
    {
        sleep(5);
    }
    print("done waiting\n");
    my $eval = `perl $scriptDir/evaluate-modified.pl $outBaseDir/$currentRun $chunkDir`;
    my $fscore = getFScore($eval);
    print("current fscore is $fscore\n");

    system("rm -r $outBaseDir/$currentRun"); #clean up after each run of Decoder

    $currentRun++;

    $dpHash{$key} = $fscore;

    return $fscore;
}

#params: param1, new-value1, param2, new-value2, etc (key-value pairs)
sub changeParams
{
    my %paramHash = ();
    while(@_ > 0)
    {
        my $inParam = shift;
        my $inVal = shift;
        $paramHash{$inParam} = $inVal;
    } 

    open(OLD_PARAM_FILE, $decoderParamFile);
    my @output = <OLD_PARAM_FILE>;
    close(OLD_PARAM_FILE);
    open(NEW_PARAM_FILE, "> $decoderParamFile");
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

#params: evaluate-output(string holding the output of the perl script)
sub getFScore
{
    my $eval = shift;
    my @eval = split(/\n/, $eval);
    my $resultLine = $eval[7];
    die "resultLine not defined" unless defined($resultLine);
    chomp($resultLine);
    my ($name, $n1, $n2, $n3, $sp1, $sp2, $st1, $st2, $e1, $e2, $e3, $g1, $g2)
        = ($resultLine =~ /(\w+)\s*\|(\d+)%\s+(\d+)%\s+(\d+)%\|(\d+)%\s+(\d+)%\|(\d+)%\s+(\d+)%\|(\d+)%\s+(\d+)%\s+(\d+)%\|(\d+)%\s+\((\d+)\)/);
    #regex must match: 'Decoder  |92% 99%  95%|75% 82%|67% 63%|70% 74% 72%|37% (154)'
    return $e3;
}
