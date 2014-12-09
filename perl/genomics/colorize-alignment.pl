#!/usr/bin/perl
use strict;
use GffTranscriptReader;

my $usage="$0 <alignment.txt> <query.gff> <subject.gff>";
die "$usage\n" unless @ARGV==3;
my ($alignmentFile,$queryGffFile,$subjectGffFile)=@ARGV;

my $reader=new GffTranscriptReader;
my $queryTranscripts=$reader->loadGFF($queryGffFile);
my $subjectTranscripts=$reader->loadGFF($subjectGffFile);
my $queryTranscript=$queryTranscripts->[0];
my $subjectTranscript=$subjectTranscripts->[0];

my $debug=@$queryTranscripts;
print "$debug query transcripts\n";

my $numQueryExons=$queryTranscript->numExons();
my $numSubjectExons=$subjectTranscript->numExons();

my $red="<span style=\"color: rgb(255, 0, 0);\">";
my $blue="<span style=\"color: rgb(51, 51, 255);\">";
my $green="<span style=\"color: rgb(51, 255, 51);\">";
my $purple="<span style=\"color: rgb(204, 51, 204);\">";
my $orange="<span style=\"color: rgb(255, 204, 0);\">";
my $pink="<span style=\"color: rgb(255, 153, 255);\">";
my $black="<span style=\"color: rgb(0, 0, 0);\">";
my $grey="<span style=\"color: rgb(192, 192, 192);\">";
my @colors=($red,$blue,$green,$purple,$orange,$pink,$black,$grey);
my $numColors=@colors;

my $queryColorIndex=0;
my $subjectColorIndex=0;
my $queryExonIndex=0;
my $subjectExonIndex=0;
my $queryExon=$queryTranscript->getIthExon(0);
my $subjectExon=$subjectTranscript->getIthExon(0);
my $queryExonLen=$queryExon->getLength()/3;
my $subjectExonLen=$subjectExon->getLength()/3;
my $queryExonConsumed=0;
my $subjectExonConsumed=0;

print "<html><body style=\"color: rgb(0, 0, 0); background-color: rgb(255, 255, 255);\"
 link=\"#000099\" vlink=\"#990099\" alink=\"#000099\">\n";
print "<pre><big>\n";
open(IN,$alignmentFile) || die "can't open file \"$alignmentFile\"\n";
while(<IN>)
  {
    if(/Query:\s*(\S+)/)
      {
	my $seq=$1;
	my $len=(length $seq);
	print "Query: ";
	my $color=$colors[$queryColorIndex%$numColors];
	print "$color";
	for(my $i=0 ; $i<$len ; ++$i)
	  {
	    my $residue=substr($seq,$i,1);
	    print "$residue";
	    if($residue ne '-' && ++$queryExonConsumed>=$queryExonLen)
	      {
		$queryExonConsumed-=$queryExonLen;
		$color=$colors[(++$queryColorIndex)%$numColors];
		$queryExon=$queryTranscript->getIthExon(++$queryExonIndex);
		$queryExonLen=($queryExon ? $queryExon->getLength()/3 : 1);
		print "</span>$color";
	      }
	  }
	print "</span>\n";
      }
    elsif(/Sbjct:\s*(\S+)/)
      {
	my $seq=$1;
	my $len=(length $seq);
	print "Sbjct: ";
	my $color=$colors[$subjectColorIndex%$numColors];
	print "$color";
	for(my $i=0 ; $i<$len ; ++$i)
	  {
	    my $residue=substr($seq,$i,1);
	    print "$residue";
	    if($residue ne '-' && ++$subjectExonConsumed>=$subjectExonLen)
	      {
		$subjectExonConsumed-=$subjectExonLen;
		$color=$colors[(++$subjectColorIndex)%$numColors];
		$subjectExon=$subjectTranscript->getIthExon(++$subjectExonIndex);
		$subjectExonLen=($subjectExon ? $subjectExon->getLength()/3 : 1);
		print "</span>$color";
	      }
	  }
	print "</span>\n";
      }
    elsif(/---/) {print "<hr>\n"}
    elsif(/alignment>/) {}
    else
      {
	print;
      }
  }
close(IN);

print "</big></pre></body></html>\n";






