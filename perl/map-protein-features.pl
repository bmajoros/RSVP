#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use GffReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <transcript.gff> <domains.gff>\n" unless @ARGV==2;
my ($transcriptFile,$domainFile)=@ARGV;

my $tReader=new GffTranscriptReader();
my $transcriptArray=$tReader->loadGFF($transcriptFile);
my $numTranscripts=@$transcriptArray;
die "wrong number of transcripts: $numTranscripts\n" unless $numTranscripts==1;
my $transcript=$transcriptArray->[0];
my $numExons=$transcript->numExons();

my $reader=new GffReader();
my $domains=$reader->loadGFF($domainFile);
foreach my $domain (@$domains) {
  my $pBegin=$domain->getBegin();
  my $pEnd=$domain->getEnd();
  my $transBegin=3*$pBegin;
  my $transEnd=3*$pEnd;
  my $dnaBegin=mapit($transBegin);
  my $dnaEnd=mapit($transEnd);
  $domain->setBegin($dnaBegin);
  $domain->setEnd($dnaEnd);
  my $gff=$domain->toGff();
  print "$gff";
}

sub mapit {
  my ($pos)=@_;
  my $upto=0;
  for(my $i=0 ; $i<$numExons ; ++$i) {
    my $exon=$transcript->getIthExon($i);
    my $strand=$exon->getStrand();
    my $exonLen=$exon->getLength();
    $upto+=$exonLen;
    if($upto>$pos) {
      $upto-=$exonLen;
      my $offset=$pos-$upto;
      my $dnaPos=$strand eq "+" ? $exon->getBegin()+$offset
	: $exon->getEnd()-1-$offset;
      return $dnaPos;
    }
  }
  die "$pos > $upto";
}





