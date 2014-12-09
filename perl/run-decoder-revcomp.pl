#!/usr/bin/perl
use strict;
use ConfigFile;
#A wrapper script for run-decoder.pl that creates two gff outputs,
#one for the forward strand and one for the reverse strand. This
#wrapper calls genezilla, gets the read profiles, and calls run-decoder.pl,
#revcomping the input and output for forward strand prediction.

#input: chunk-prefix, out-prefix, param-file

die "run-decoder-revcomp.pl <data-prefix> <out-prefix> <param-file>" unless @ARGV == 3;

my $dataPrefix = $ARGV[0];
my $outPrefix = $ARGV[1];
my $paramFile = $ARGV[2];

my $configFile=new ConfigFile($paramFile);

my $packageDir=$configFile->lookup("package");
my $genezillaDir = $configFile->lookup("genezilladir");
my $iso = $configFile->lookup("iso");
my $scriptDir = $configFile->lookup("scriptdir");
my $bamFile = $configFile->lookup("bamfile");
my $samtoolsPath = $configFile->lookup("samtoolsdir");
my $billScriptDir = $configFile->lookup("billscriptdir");
my $strandedReads = $configFile->lookup("strandedreads");
my $runBothStrands = $configFile->lookup("runbothstrands");
my $usePileup = $configFile->lookup("genezilla_use_pileup");
my $chunksDir = $configFile->lookup("chunkdir");
my $dashP=$usePileup ? "-P $chunksDir" : "";

#get coordinate and strand info
open(COORD, "$dataPrefix.coords.gff");
my $str = <COORD>;
my @split = split(/\s+/,$str);
my $chr = $split[0];
my $begin = $split[3];
my $end = $split[4];
my $strand = $split[6]; #if we know the strand we want to look at, this will be + or -

die "Error: either run on both strands or specify transcript strand in coords file" if $strand eq "." && $runBothStrands eq "false";

#reverse strand
if ($strand eq "-" || $runBothStrands eq "true")
{
    #get genezilla graph
    system("$genezillaDir/genezilla $dashP $iso $dataPrefix.fasta -o $dataPrefix.graph > $dataPrefix.genezilla; gzip -f --best $dataPrefix.graph");
    die unless -e "$dataPrefix.genezilla";
    die unless -e "$dataPrefix.graph.gz";

    #get read profiles
    if ($strandedReads eq "true")
    {
        system("perl $scriptDir/get-read-profiles.pl $chr $begin $end - - $dataPrefix $bamFile $samtoolsPath");
    }
    else
    {
        system("perl $scriptDir/get-read-profiles.pl $chr $begin $end 0 - $dataPrefix $bamFile $samtoolsPath");
    }

    #call run-decoder.pl
    system("perl $scriptDir/run-decoder.pl $dataPrefix $outPrefix.reverse $paramFile");
}

#forward strand
if ($strand eq "+" || $runBothStrands eq "true")
{
    #get genezilla graph
    system("perl $billScriptDir/revcomp-fasta.pl  $dataPrefix.fasta > $dataPrefix.revcomp.fasta");
    system("$genezillaDir/genezilla $dashP $iso $dataPrefix.revcomp.fasta -o $dataPrefix.revcomp.graph > $dataPrefix.revcomp.genezilla; gzip -f --best $dataPrefix.revcomp.graph");

    #get forward strand, revcomped read profiles
    if ($strandedReads eq "true")
    {
        system("perl $scriptDir/get-read-profiles.pl $chr $begin $end + + $dataPrefix.forward $bamFile $samtoolsPath");
    }
    else
    {
        system("perl $scriptDir/get-read-profiles.pl $chr $begin $end 0 + $dataPrefix.forward $bamFile $samtoolsPath");
    }   


    system("perl $scriptDir/revcomp-profile.pl $dataPrefix.forward.readProfile > $dataPrefix.revcomp.readProfile");
    system("perl $scriptDir/revcomp-profile.pl $dataPrefix.forward.spliceInProfile > $dataPrefix.revcomp.spliceInProfile");
    system("perl $scriptDir/revcomp-profile.pl $dataPrefix.forward.spliceOutProfile > $dataPrefix.revcomp.spliceOutProfile");

    #get revcomped splice-read profile
    system("perl $scriptDir/revcomp-splice-profile.pl $dataPrefix.spliceProfile > $dataPrefix.revcomp.spliceProfile");

    #call run-decoder.pl
    system("perl $scriptDir/run-decoder.pl $dataPrefix.revcomp $outPrefix.revcomp.reverse $paramFile");

    #revcomp output file
    system("perl $scriptDir/revcomp-gff.pl $outPrefix.revcomp.reverse.gff > $outPrefix.forward.gff");
}


#if running on both strands, merge predictions into one gff file
if ($runBothStrands eq "true")
{
    system("perl $scriptDir/merge-predictions.pl $outPrefix $dataPrefix > $outPrefix.gff");
}
elsif($strand eq "-")
{
    system("cp $outPrefix.reverse.gff $outPrefix.gff");
}
elsif($strand eq "+")
{
    system("cp $outPrefix.forward.gff $outPrefix.gff");
}

#remove parses that scored -Infinity
system("perl $scriptDir/remove-inf.pl $outPrefix.gff > $outPrefix.temp.gff");
system("mv $outPrefix.temp.gff $outPrefix.gff");

#cleanup
#exit;###
system("rm $dataPrefix.revcomp.fasta");
system("rm $dataPrefix.genezilla");
#system("rm $dataPrefix.graph.gz");
system("rm $dataPrefix.revcomp.genezilla");
system("rm $dataPrefix.revcomp.graph.gz");
system("rm $dataPrefix.readProfile");
system("rm $dataPrefix.spliceInProfile");
system("rm $dataPrefix.spliceOutProfile");
system("rm $dataPrefix.revcomp.readProfile");
system("rm $dataPrefix.revcomp.spliceInProfile");
system("rm $dataPrefix.revcomp.spliceOutProfile");
system("rm $dataPrefix.revcomp.spliceProfile");
system("rm $dataPrefix.forward.readProfile");
system("rm $dataPrefix.forward.spliceInProfile");
system("rm $dataPrefix.forward.spliceOutProfile");
system("rm $outPrefix.revcomp.reverse.gff");
system("rm $outPrefix.reverse.gff");
system("rm $outPrefix.forward.gff");
