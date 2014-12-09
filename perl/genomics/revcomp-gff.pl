#!/usr/bin/perl
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffReader;
$|=1;

#########################################################################
# Computes the reverse-complement of a GFF file and outputs the
# result as a new GFF file on STDOUT.
#
# bmajoros@tigr.org
#########################################################################

my $usage="$0 <*.gff> <substrate-length>\n";
die $usage unless @ARGV==2 || @ARGV==1;
my ($filename,$seqLen)=@ARGV;

my $reader=new GffReader;
my $features=$reader->loadGFF($filename);

my $n=@$features;
if(!defined($seqLen))
  {
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $feature=$features->[$i];
	my $begin=$feature->getBegin();
	my $end=$feature->getEnd();
	if($end>$seqLen) {$seqLen=$end}
      }
  }
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $feature=$features->[$i];
    my $begin=$feature->getBegin();
    my $end=$feature->getEnd();
    $feature->setBegin($seqLen-$end);
    $feature->setEnd($seqLen-$begin);
    $feature->{strand}=comp($feature->{strand});
    my $output=$feature->toGff();
    print $output;
  }

#------------------------------------------------------
sub comp
  {
    my ($strand)=@_;
    if($strand eq "+") {return "-"}
    if($strand eq "-") {return "+"}
    if($strand eq ".") {return "."}
    die "Unknown strand \"$strand\" in $0";
  }





