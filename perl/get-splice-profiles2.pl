#!/usr/bin/perl

use warnings;
use strict;

die "get-splice-profiles.pl <chunk-dir> <intron-file>" unless @ARGV==2;
my ($chunkDir, $bedFile) = @ARGV;

open(BED, $bedFile);
my %spliceHash;
while (<BED>) {
  my @lineSplit=split(/\s+/, $_);
  my ($chr, $begin, $end, $num) = 
    ($lineSplit[0], $lineSplit[1], $lineSplit[2], $lineSplit[3]);
  push @{$spliceHash{$chr}},[$begin,$end,$num];
}
close(BED);


my @keys=keys %spliceHash;
foreach my $key (@keys) {
  $array=$spliceHash{$key};
  @array=sort {$a->[0] <=> $b->[0]} @$array;
}

my @chunkList = `ls $chunkDir/*.coords.gff`;
my $numChunks=@chunkList;
for(my $i=0 ; $i<$numChunks ; ++$i) {
  my $coordFile=$chunkList[$i];
  my ($chunkNum) = ($coordFile =~ /(\d+)\.coords\.gff/);
  my $spliceProfile = "$chunkDir/$chunkNum.spliceProfile";
  open(COORDS, $coordFile);
  my $line = <COORDS>;
  close(COORDS);
  chomp($line);
  my @lineSplit = split(/\s+/, $line);
  my ($chr, $begin, $end) = ($lineSplit[0], $lineSplit[3], $lineSplit[4]);
  my $header = "@ $chr $begin $end";
  open(OUT, ">$spliceProfile");
  print(OUT "$header\n");
  my $records=$spliceHash{$chr};
  my $N=@$records;
  my $from=0;
  my $to=$N;
  while($from<$to) {
    my $mid=($from+$to)/2;
    my $rec=$records->[$mid];
    my ($donor,$acceptor,$count)=@$rec;
    if($end<$donor) { $to=$mid-1 }
    else if($begin>$acceptor) { $from=$mid+1 }
    else { last }
  }
  my $j=$from;
  for( ; $j>0 ; --$j) {
    if($records->[$j]->[1]<$begin) { last }
  }
  for( ; $j<numRecords ; ++$j) {
    if($records->[$j]->[0]>$end) { last }
    my ($donor,$acceptor,$count)=@{$records->[$j]};
    if($donor>$begin && $acceptor<$end) {
      my $relBegin=$donor-$begin;
      my $relEnd=$acceptor-$begin;
      print OUT "$relBegin#$relEnd\t";
    }
  }

===
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
  close(OUT);
}
