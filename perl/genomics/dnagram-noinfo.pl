#!/usr/bin/perl
use strict;
use lib('/home/ohler/bmajoros/genomics/browser',
        '/home/ohler/bmajoros/genomics/perl',
        '/home/ohler/bmajoros/perlib',
        '.');
use GeneScript::Color;
use GeneScript::Signal;
use GeneScript::Canvas;
use GeneScript::PostscriptPage;
use GeneScript::ColorLegend;
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile> <outfile>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my ($numSeqs,@sampleSizes);
my $reader=new FastaReader($infile);
my (@counts);
while(1)
  {
    my ($def,$seq)=$reader->nextSequence();
    last unless $def;
    ++$numSeqs;
    my $len=length($seq);
    for(my $i=0 ; $i<$len ; ++$i)
      {
	my $letter=substr($seq,$i,1);
	next if($letter eq "-");
	++$counts[$i]->{$letter};
	++$sampleSizes[$i];
      }
  }
my @information;
my $len=@counts;
for(my $i=0 ; $i<$len ; ++$i)
  {
    my $hash=$counts[$i];
    my @keys=keys %$hash;
    my $sum=0;
    foreach my $letter (@keys) {$sum+=$hash->{$letter}}
    my $H=0;
    foreach my $letter (@keys) {
      #print "$letter $hash->{$letter} ";
      $hash->{$letter}/=$sum;
      my $p=$hash->{$letter};
      $H-=$p*log($p)/log(2) unless $p==0.0;
    }
    #print "\n";
    my $maxH=log(4)/log(2);
    my $I=($maxH-$H)/$maxH;
    $information[$i]=1;
  }

my @colors=($Color::red,$Color::blue,$Color::green,$Color::gray,
	    $Color::brown,$Color::purple,$Color::lightGreen,
	    $Color::lightBlue,$Color::hotPink,$Color::gold,
	    $Color::darkGray,$Color::yellow,$Color::orange);
my %colorMap=%{{"A"=>Color::fromGrayscale(0.0),   #$Color::red,
	   "T"=>Color::fromGrayscale(0.3),        #$Color::purple,
	   "C"=>Color::fromGrayscale(0.55),        #$Color::green,
           "G"=>Color::fromGrayscale(0.85),        #$Color::orange}};
	    }};
my %colorMap=%{{"A"=>$Color::red,
	   "T"=>$Color::purple,
	   "C"=>$Color::green,
           "G"=>$Color::orange
	    }};

my $page=new GeneScript::PostscriptPage;
my $canvas=new GeneScript::Canvas;

my $left=0.1;
my $right=0.9;
my $xRange=$right-$left;
my $colWidth=$xRange/$len;
my $x=$left;
my $bottom=0.7;
my $top=0.9;
my $yRange=$top-$bottom;
for(my $i=0 ; $i<$len ; ++$i)
  {
    my $hash=$counts[$i];
    my @keys=keys %$hash;
    @keys=sort {$hash->{$b} <=> $hash->{$a}} @keys;
    my $y=$top;
    $y-=(1-$information[$i])*$yRange;
    #$y-=(1-$information[$i]*$sampleSizes[$i]/$numSeqs)*$yRange;
    foreach my $letter (@keys)
      {
	my $color=$colorMap{$letter};
	if(!defined($color)) {$color=$Color::black}
	$canvas->setColor($color);
	my $P=$hash->{$letter};
	my $height=$P*$yRange;
	my $fudgeFactor;
	if($letter eq "A") {$fudgeFactor=1.36}
	if($letter eq "C") {$fudgeFactor=1.36}
	if($letter eq "G") {$fudgeFactor=1.30}
	if($letter eq "T") {$fudgeFactor=1.36}
	my $yAdj=0;
	if($letter eq "G") {$yAdj=0.03*$height} #0.003+($height-0.1)*0.03}
	$height*=$information[$i];
	#my $sampleSizeRatio=$sampleSizes[$i]/$numSeqs;
	#$height*=$sampleSizeRatio;
	$canvas->drawText($letter,$x,$y-$height+$yAdj,$colWidth,
			  $fudgeFactor*$height);
	$y-=$height;
      }
    $x+=$colWidth;
  }
$canvas->setColor($Color::black);
$canvas->setLineWidth(0.001);
$canvas->drawRect($left,$bottom,$xRange,$yRange);

$page->attach($canvas,0,0,1,1);
$page->changeMargins(0,0);
my $ps=$page->toPostscript();
open(PS,">$outfile") || die;
print PS $ps;
close(PS);

# Print WMM onto STDOUT
foreach my $letter ("A","C","G","T")
  {
    print "$letter\t";
    for(my $i=0 ; $i<$len ; ++$i)
      {
	my $hash=$counts[$i];
	my $P=$hash->{$letter};
	$P=int(1000*$P+5/9)/1000;
	print "$P\t";
      }
    print "\n";
  }







