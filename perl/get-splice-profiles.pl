#!/usr/bin/perl

use warnings;
use strict;

die "get-splice-profiles.pl <chunk-dir> <intron-file>" unless @ARGV==2;

my ($chunkDir, $bedFile) = @ARGV;

open(BED, $bedFile);
my %spliceHash = ();
while (<BED>)
{
  my @lineSplit = split(/\s+/, $_);
  my ($chr, $begin, $end, $num) = ($lineSplit[0], $lineSplit[1], $lineSplit[2], $lineSplit[3]);
  my $key = "$chr:$begin#$end";
  $spliceHash{$key} = $num;
}
close(BED);

my $chunks = `ls $chunkDir/*.coords.gff`;

my @chunkList = split(/\s+/, $chunks);
for my $coordFile (@chunkList)
{
  my ($chunkNum) = ($coordFile =~ /(\d+)\.coords\.gff/);
  my $spliceProfile = "$chunkDir/$chunkNum.spliceProfile";
  open(COORDS, $coordFile);
  open(OUT, "> $spliceProfile");
  my $line = <COORDS>;
  chomp($line);
  my @lineSplit = split(/\s+/, $line);
  my ($chr, $begin, $end) = ($lineSplit[0], $lineSplit[3], $lineSplit[4]);
  my $header = "@ $chr $begin $end";
  print(OUT "$header\n");
  for my $key (keys(%spliceHash))
  {
    my ($keyChr, $keyBegin, $keyEnd) = ($key =~ /(\w+):(\d+)#(\d+)/);
    if ($keyBegin >= $begin && $keyEnd <= $end && $keyChr eq $chr)
    {
      my $relBegin = $keyBegin - $begin; #using 0-based coordinates
      my $relEnd = $keyEnd - $begin;
      my $newKey = "$relBegin#$relEnd";
      print(OUT "$newKey\t$spliceHash{$key}\n");
    }
  }
  close(COORDS);
  close(OUT);
}
