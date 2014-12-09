#!/usr/bin/perl
use strict;
use GffReader;
use GffTranscriptReader;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <hypothetical-gff> <correct-gff>";
die "$usage\n" unless @ARGV==2;
my ($filename1,$filename2)=@ARGV;

my $reader=new GffReader;
my $hypothFeatures=$reader->loadGFF($filename1);
my $correctFeatures=$reader->loadGFF($filename2);

my $leftTerminus=getLeftTerminus($correctFeatures);
my $rightTerminus=getRightTerminus($correctFeatures);
my $length=$rightTerminus-$leftTerminus;

my $hypothVector=makeVector($hypothFeatures);
my $correctVector=makeVector($correctFeatures);
my ($percentIdentity,$numDifferences)=
  getPercentIdentity($hypothVector,$correctVector);
my ($sens,$spec,$TP,$TN,$FP,$FN)=
  getSensAndSpecificity($hypothVector,$correctVector);
$percentIdentity=int(1000*$percentIdentity)/10;
$sens=int(1000*$sens)/10;
$spec=int(1000*$spec)/10;
print "$percentIdentity\% identity ($numDifferences differences out of $length), sensitivity=$sens\%, specificity=$spec\%, TP=$TP TN=$TN FP=$FP FN=$FN\n";
my $firstDifference=findFirstDifference($hypothVector,$correctVector);
if(defined($firstDifference)) {print "First difference at $firstDifference\n"}
#-----------------------------------------------------------------
sub findFirstDifference
  {
    my ($vector1,$vector2)=@_;
    #my $n=max(0+@$vector1,0+@$vector2);
    for(my $i=$leftTerminus ; $i<$rightTerminus ; ++$i)
      {
	my $elem1=$vector1->[$i];
	my $elem2=$vector2->[$i];
	if($elem1!=$elem2){return $i}
      }
    return undef;
  }
#-----------------------------------------------------------------
sub max
  {
    my ($a,$b)=@_;
    return ($a>$b ? $a : $b);
  }
#-----------------------------------------------------------------
sub getPercentIdentity
  {
    my ($vector1,$vector2)=@_;
    my $identity;
    for(my $i=$leftTerminus ; $i<$rightTerminus ; ++$i)
      {
	my $elem1=$vector1->[$i];
	my $elem2=$vector2->[$i];
	if($elem1==$elem2){++$identity}
      }
    return ($identity/$length,$length-$identity);
  }
#-----------------------------------------------------------------
sub getSensAndSpecificity
  {
    my ($hypoth,$correct)=@_;
    my ($TP,$TN,$FP,$FN)=(0,0,0,0);
    for(my $i=$leftTerminus ; $i<$rightTerminus ; ++$i)
      {
	my $hypothElem=$hypoth->[$i];
	my $correctElem=$correct->[$i];
	if($hypothElem==$correctElem)
	  {
	    if($correctElem==1){++$TP}
	    else {++$TN}
	  }
	else
	  {
	    if($correctElem==1){++$FN}
	    else {++$FP}
	  }
      }
    my $sensitivity=($TP+$FN>0 ? $TP/($TP+$FN) : 1);
    my $specificity=($TP+$FP>0 ? $TP/($TP+$FP) : 1);
    return ($sensitivity,$specificity,$TP,$TN,$FP,$FN);
  }
#-----------------------------------------------------------------
sub makeVector
  {
    my ($features)=@_;
    my $vector=[];
    my $n=@$features;
    for(my $i=0 ; $i<$n ; ++$i) {$vector->[$i]=0}
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $feature=$features->[$i];
	next unless $feature->{featureType}=~/exon/;
	my $begin=$feature->getBegin()-1;
	my $end=$feature->getEnd();
	if($begin<0) { $begin=0 }
	for(my $j=$begin ; $j<$end ; ++$j)
	  {
	    $vector->[$j]=1;
	  }
      }
    return $vector;
  }
#-----------------------------------------------------------------
#my $leftTerminus=getLeftTerminus($correctFeatures);
sub getLeftTerminus
  {
    my ($features)=@_;
    my $left;
    foreach my $feature (@$features)
      {
	if(!defined($left) || $feature->getBegin()<$left)
	  { $left=$feature->getBegin() }
      }
    return $left;
  }
#-----------------------------------------------------------------
#my $rightTerminus=getRightTerminus($correctFeatures);
sub getRightTerminus
  {
    my ($features)=@_;
    my $right;
    foreach my $feature (@$features)
      {
	if(!defined($right) || $feature->getEnd()>$right)
	  { $right=$feature->getEnd() }
      }
    return $right;
  }
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#-----------------------------------------------------------------








