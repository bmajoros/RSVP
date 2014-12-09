#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffReader;
use FastaReader;
use CodonIterator;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename,$fastaFilename)=@ARGV;

my $reader=new GffReader();
my $featureArray=$reader->loadGFF($filename);
my $numFeatures=@$featureArray;
my $substrate=($numFeatures>0 ? $featureArray->[0]->{substrate} : ".");

for(my $i=0 ; $i<$numFeatures ; ++$i)
  {
    my $feature=$featureArray->[$i];
    my $strand=$feature->{strand};
    my $type=$feature->getType();
    $type="\L$type";
    if($strand eq "+")
      {
	if($type eq "single-exon")
	  {
	    #$feature->setBegin($feature->getBegin()-3);
	    #$feature->setEnd($feature->getEnd()+3);
	  }
	elsif($type eq "initial-exon")
	  {
	    #$feature->setBegin($feature->getBegin()-3);
	  }
	elsif($type eq "internal-exon")
	  {
	  }
	elsif($type eq "final-exon")
	  {
	    #$feature->setEnd($feature->getEnd()+3);
	  }
	else {next}
      }
    else
      {
	if($type eq "single-exon")
	  {
	    #$feature->setBegin($feature->getBegin()-3);
	    #$feature->setEnd($feature->getEnd()+3);
	  }
	elsif($type eq "initial-exon")
	  {
	    #$feature->setEnd($feature->getEnd()+3);
	  }
	elsif($type eq "internal-exon")
	  {
	  }
	elsif($type eq "final-exon")
	  {
	    #$feature->setBegin($feature->getBegin()-3);
	  }
	else {next}
      }
    my $score=$feature->getScore();
    my $extra=$feature->{additionalFields};
    if($score=~/,/)
      {
	my @scores=split/,/,$score;
	$score="";
	my $numScores=0;
	foreach my $s (@scores)
	  {if($s=~/\d/){$score.=(length($score)?",$s":$s);++$numScores}}
	$feature->setScore($score);
      }
    my $gff=$feature->toGff();
    print "$gff";
  }


