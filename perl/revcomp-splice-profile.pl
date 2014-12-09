#!/usr/bin/perl
use warnings;
use strict;

die "revcomp-splice-profile.pl <*.spliceProfile>" unless @ARGV==1;

my $file = shift;
open(SP, $file);

my @sp = <SP>;
my $header = $sp[0];
print("$header");
my @headerSplit = split(/\s+/, $header);
my ($chunkBegin, $chunkEnd) = ($headerSplit[2], $headerSplit[3]);
my $length = $chunkEnd - $chunkBegin + 1;
for (my $i = 1; $i < @sp; $i++)
{
  my @lineSplit = split(/\s+/, $sp[$i]);
  my $oldKey = $lineSplit[0];
  my $val = $lineSplit[1];
  my ($oldBegin, $oldEnd);
  if ($oldKey =~ /(\d+)#(\d+)/)
  {
    ($oldBegin, $oldEnd) = ($1, $2);
  }
  else
  {
    die "error: key formatted incorrectly";
  }
  my $newBegin = $length - $oldEnd - 1;
  my $newEnd = $length - $oldBegin - 1;
  my $newKey = "$newBegin#$newEnd";
  print("$newKey\t$val\n");
}
