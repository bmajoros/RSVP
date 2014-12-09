#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffTranscriptReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffTranscriptReader();
my $transcriptArray=$reader->loadGFF($filename);
my $numTranscripts=@$transcriptArray;
my ($substrate,$source,$strand);
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcriptArray->[$i];
    $strand=$transcript->{strand};
    my $numExons=$transcript->numExons();
    $substrate=$transcript->{substrate};
    $source=$transcript->{source};

    if($strand eq "+")
      {
	#####################
	# HANDLE INITIAL EXON
	#####################
	my $initialExon=$transcript->getIthExon(0);
	printATG($initialExon);
	printGT($initialExon);
	
	#######################
	# HANDLE INTERNAL EXONS
	#######################
	for(my $j=1 ; $j<$numExons-1 ; ++$j)
	  {
	    my $exon=$transcript->getIthExon($j);
	    printAG($exon);
	    printGT($exon);
	  }
	
	###################
	# HANDLE FINAL EXON
	###################
	my $finalExon=$transcript->getIthExon($numExons-1);
	printAG($finalExon);
	printTAG($finalExon);
      }
    else
      {
	###################
	# HANDLE FINAL EXON
	###################
	my $finalExon=$transcript->getIthExon($numExons-1);
	printTAG($finalExon);
	printAG($finalExon);

	#######################
	# HANDLE INTERNAL EXONS
	#######################
	for(my $j=$numExons-2 ; $j>1 ; --$j)
	  {
	    my $exon=$transcript->getIthExon($j);
	    printGT($exon);
	    printAG($exon);
	  }
	
	#####################
	# HANDLE INITIAL EXON
	#####################
	my $initialExon=$transcript->getIthExon(0);
	printGT($initialExon);
	printATG($initialExon);
      }
  }


sub printATG
  {
    my ($initialExon)=@_;
    my ($atgBegin, $atgEnd);
    if($strand eq '+') {$atgBegin=$initialExon->{begin}+1}
    else {$atgBegin=$initialExon->{end}-2}
    $atgEnd=$atgBegin+2;
    print "$substrate\tvertex\tATG\t$atgBegin\t$atgEnd\t.\t$strand\t.\n";
  }

sub printTAG
  {
    my ($finalExon)=@_;
    my ($tagBegin, $tagEnd);
    if($strand eq "+") {$tagBegin=$finalExon->{end}-2}
    else {$tagBegin=$finalExon->{begin}+1}
    $tagEnd=$tagBegin+2;
    print "$substrate\tvertex\tTAG\t$tagBegin\t$tagEnd\t.\t$strand\t.\n";
  }

sub printGT
  {
    my ($exon)=@_;
    my ($gtBegin, $gtEnd);
    if($strand eq "+")
      {
	$gtBegin=$exon->{end}+1;
	$gtEnd=$gtBegin+1;
	print "$substrate\tvertex\tGT\t$gtBegin\t$gtEnd\t.\t$strand\t.\n";
      }
    else
      {
	$gtBegin=$exon->{begin}-1;
	$gtEnd=$gtBegin+1;
	print "$substrate\tvertex\tGT\t$gtBegin\t$gtEnd\t.\t$strand\t.\n";
      }
  }

sub printAG
  {
    my ($exon)=@_;
    my ($agBegin, $agEnd);
    if($strand eq "+")
      {
	$agBegin=$exon->{begin}-1;
	$agEnd=$agBegin+1;
	print "$substrate\tvertex\tAG\t$agBegin\t$agEnd\t.\t$strand\t.\n";
      }
    else
      {
	$agBegin=$exon->{end}+1;
	$agEnd=$agBegin+1;
	print "$substrate\tvertex\tAG\t$agBegin\t$agEnd\t.\t$strand\t.\n";
      }
  }

