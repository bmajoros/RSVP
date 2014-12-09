#!/usr/bin/perl
use ConfigFile;
use strict;

#Wrapper script for run-decoder-revcomp.pl that uses the parameter file for paths.
#Builds a .q wrapper file for run-decoder-revcomp.pl for each chunk

die "usage: batch-decoder-revcomp.pl <param-file>"  unless @ARGV == 1;

my $paramFile = $ARGV[0];
my $maxJobs = 1000;

my $configFile=new ConfigFile($paramFile);

my $dump=$configFile->lookup("dump");
die "error: do not run batch jobs with dump on" if $dump eq "true";

my $chunkDir=$configFile->lookup("chunkdir");
my $outDir=$configFile->lookup("outdir");
my $scriptDir=$configFile->lookup("scriptdir");
my $usePileup = $configFile->lookup("genezilla_use_pileup");
my $pileupDir = $configFile->lookup("pileups");
my $packageDir=$configFile->lookup("package");

# Generate pileup files for individual chunks for GeneZilla
if($usePileup) {
  my $cmd="$packageDir/perl/get-chunk-pileups.pl $chunkDir $pileupDir $packageDir";
  system($cmd);
}

my $files = `ls $chunkDir/*.fasta`;
my @files = split(/\s+/, $files);
my @prefixes = ();
for my $file (@files)
{
    my ($filePrefix) = ($file =~ /(\w+).fasta/);
    push(@prefixes, $filePrefix);
}

while (@prefixes > 0)
{
    my $numJobs = `qstat | wc -l`;  #throttling jobs
    while ($numJobs >= $maxJobs)
    {
        my @jobs = `qstat`; #delete 'Eqw' jobs and put the chunk back in the list
        for (my $j = 2; $j < @jobs; $j++) #two lines of header info in qstat output
        {
            my $job = $jobs[$j];
            my @jobSplit = split(/\s+/, $job);
            my $jobID = $jobSplit[0];
            my $jobName = $jobSplit[2];
            my ($chunkID) = ($jobName =~ /c(\w+)\.q/);
            my $jobStatus = $jobSplit[4];
            if ($jobStatus eq "Eqw")
            {
                system("qdel $jobID");
                push(@prefixes, $chunkID);
            }
        }
        sleep(30);
        $numJobs = `qstat | wc -l`;
    }

    my $prefix = shift(@prefixes);
    my $chunkPrefix = "$chunkDir/$prefix";
    my $outPrefix = "$outDir/$prefix";
    my $qfile = makeQsub($chunkPrefix, $outPrefix, $paramFile, $scriptDir);
    system("qsub $qfile");
    system("rm $qfile");
}

my $numJobs = `qstat | wc -l`;
while ($numJobs > 0) #mostly a duplicate of above code
{
    my @jobs = `qstat`;
    for (my $j = 2; $j < @jobs; $j++) #two lines of header info in qstat output
    {
        my $job = $jobs[$j];
        my @jobSplit = split(/\s+/, $job);
        my $jobID = $jobSplit[0];
        my $jobName = $jobSplit[2];
        my ($chunkID) = ($jobName =~ /c(\w+)\.q/);
        my $jobStatus = $jobSplit[4];
        if ($jobStatus eq "Eqw") #different from above here
        {
            system("qdel $jobID");
            my $ecPrefix = "$chunkDir/$chunkID";
            my $eoPrefix = "$outDir/$chunkID";
            my $eqfile = makeQsub($ecPrefix, $eoPrefix, $paramFile, $scriptDir);
            system("qsub $eqfile");
            system("rm $eqfile");
        }
    }
    sleep(30);
    $numJobs = `qstat | wc -l`;
}

#system("cat $outDir/output-*.txt > $outDir/output.txt");
system("rm $outDir/output-*.txt");
system("cat $outDir/error-*.txt > $outDir/error.txt");
system("rm $outDir/error-*.txt");

#my @jobs = `qstat`;
#for (my $j = 2; $j < @jobs; $j++) #two lines of header info in qstat output
#{
#    my $job = $jobs[$j];
#    @jobSplit = split(/\s+/, $job);
#    my $jobID = $jobSplit[0];
#    my $jobName = $jobSplit[2];
#    my ($chunkID) = ($jobName =~ /c(\w+)\.q/);
#    my $jobStatus = $jobSplit[4];
#    if ($jobStatus eq "Eqw")
#    {
#        system("qdel $jobID");
#        push(@prefixes, $chunkID);
#    }
#}


sub makeQsub($chunkPrefix, $outPrefix, $paramFile, $scriptDir)
{
    my ($chunkPrefix, $outPrefix, $paramFile, $scriptDir) = @_;
    my ($chunkNum) = ($chunkPrefix =~ /(\w+)$/);
    my $outputFile = "c$chunkNum.q";
    open(OUT, "> $outputFile");
    print OUT "\#\$ -S /bin/bash -cwd\n";
    print OUT "\#\$ -o $outDir/output-$chunkNum.txt -e $outDir/error-$chunkNum.txt\n";
    print OUT "\#\$ -l highprio\n";
    print OUT "\#\$ -V\n";
    print OUT "\#\$ -N $outputFile\n";
    print OUT "perl $scriptDir/run-decoder-revcomp.pl $chunkPrefix $outPrefix $paramFile\n";
    close(OUT);
    return $outputFile;
}

