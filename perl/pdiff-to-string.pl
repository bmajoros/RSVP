#!/usr/bin/perl
#Make a pdiff file into a SAM-type string according to the format specified by Song

use warnings;
use strict;

die "pdiff-to-string.pl <pdiff-file>" unless @ARGV==1;

my $prefix = shift;

open(TEST, $prefix);

my @test = <TEST>;

my $header = $test[0];

my ($annID, $predID, $strand, $chunk) = ($header =~ /annotation: (\S+)\tprediction: (\S+)\tstrand: (\S+)\tchunk: (\S+)/)
  or die "error: pdiff header formatted incorrectly";

my ($chrNum) = ($annID =~ /AT(\d+)G/) or die "error: annotation ID does not have chromosome number or is formatted incorrectly";
my $chr = "Chr$chrNum";

my $outputString = "";
my $state;
my $stateNum;

my $firstLine = $test[1];
my @firstLineSplit = split(/\s+/, $firstLine);
my $leftMostCoord = $firstLineSplit[0];

for (my $i = 1; $i < @test; $i++)
{
  my $line = $test[$i];
  chomp ($line);
  my @lineSplit = split(/\s+/, $line);
  my ($coord, $currentState) = @lineSplit;
  if ($coord < $leftMostCoord) { die "error: found coord less than first coord"; }
  if (!(defined($state)))
  {
    $state = $currentState;
    $stateNum = 0;
  }
  if (!($currentState eq $state))
  {
    if ($outputString eq "")
    {
      $outputString = "${stateNum}${state}";
    }
    else
    {
      $outputString = "$outputString:${stateNum}${state}";
    }
    $state = $currentState;
    $stateNum = 1;
  }
  else #$currentState eq $state
  {
    $stateNum++;
  }
}

$outputString = "$outputString:${stateNum}${state}"; #add final state to string

print("$chr\t$leftMostCoord\t$strand\t$outputString\t$annID\t$predID\t$chunk\n");
