#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl');
use Exon;
use Transcript;

my $usage="\n$0 <source> <*.eff>\n(EFF=Ela's Feature Format, GFF=General Feature Format)\n\nexample:  $0 human chr21.eff\n";
die "$usage\n" unless @ARGV==2;
my ($source,$infile)=@ARGV;

my $transcriptId=0;
my $transcript;
my $substrate;
my @exons;
open(IN,$infile) || die;
while(<IN>)
  {
    if(/Sequence name: (\S+)/) {$substrate=$1}
    if(/^\s*\d+\s+\d+\s+(\S)\s+\S+\s+(\d+)\s+(\d+)/)
      {
	my ($strand,$begin,$end)=($1,$2,$3);
	my $strand="+";
	if($begin>$end) {($begin,$end)=($end,$begin); $strand="-"}
	if(!defined($transcript))
	  {
	    ++$transcriptId;
	    $transcript=new Transcript($transcriptId,$strand);
	  }
	my $exon=new Exon($begin-1,$end,$transcript);
	$transcript->addExon($exon);
	$transcript->setSubstrate($substrate);
	$transcript->setSource($source);
      }
    elsif(/^\s*$/)
      {
	next unless defined $transcript;
	$transcript->setExonTypes();
	my $gff=$transcript->toGff();
	print $gff;
	undef $transcript;
      }
  }
close(IN);

