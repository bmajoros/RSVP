#!/usr/bin/perl
use strict;
use ProgramName;
use GffTranscriptReader;

my $name=ProgramName::get();
die "$name <in.gff> <in.pileup> <out.gff> <package-dir>\n" unless @ARGV==4;
my ($inGffFile,$pileupFile,$outGffFile,$package)=@ARGV;

open(OUT,">$outGffFile") || die "Can't write to file: $outGffFile\n";
my $reader=new GffTranscriptReader();
my $transcripts=$reader->loadGFF($inGffFile);
my @pileup=`$package/dump-pileup $pileupFile`;
my $numTranscripts=@$transcripts;
my $nextID=$numTranscripts;
foreach my $transcript (@$transcripts) {
  my $gff=$transcript->toGff();
  print OUT $gff;
  my $strand=$transcript->getStrand();
  $transcript->sortExons();
  if($strand eq "-") { $transcript->reverseExonOrder() }
  my $numExons=$transcript->numExons();
  for(my $i=1 ; $i<$numExons ; ++$i) {
    my $thisExon=$transcript->getIthExon($i-1);
    my $nextExon=$transcript->getIthExon($i);
    my $intron=[$thisExon->getEnd(),$nextExon->getBegin()];
    if(covered($intron)) {
      my $newTranscript=$transcript->copy();
      $newTranscript->sortExons();
      if($strand eq "-") { $newTranscript->reverseExonOrder() }
      my $thisExon=$newTranscript->getIthExon($i-1);
      my $nextExon=$newTranscript->getIthExon($i);
      #print $thisExon->toGff();
      #print $nextExon->toGff();
      #print "================================\n";
      $thisExon->setEnd($nextExon->getEnd());
      $thisExon->setType($nextExon->getType());
      $newTranscript->deleteExon($i);
      $newTranscript->setTranscriptId($nextID++);
      my $gff=$newTranscript->toGff();
      print OUT $gff;
    }
  }
}
close(OUT);


sub covered {
  my ($intron)=@_;
  my ($begin,$end)=@$intron;
  for(my $i=$begin ; $i<$end ; ++$i) {
    if($pileup[$i]==0) { return 0 }
  }
  return 1;
}
