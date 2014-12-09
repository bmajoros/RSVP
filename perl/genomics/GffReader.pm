package GffReader;
use strict;
use Exon;
use Feature;

######################################################################
#
# GffReader.pm 
# bmajoros@tigr.org 7/8/2002
#
# Returns a list of Features, sorted by their begin coordinates.
#
# Attributes:
#
# Methods:
#   $reader=new GffReader();
#   $featureArray=$reader->loadGFF($filename);
#   $bySubstrateHash=$reader->hashBySubstrate($filename); 
#       returns a hash in which each key is a substrate ID and the 
#       corresponding value is a pointer to an array of features which 
#       are sorted by increasing coordinate
# Static methods:
#   $groupArray=groupByOverlaps($featureArray); # each group=array of features
#                                               # ^ substrate is ignored
#   $groupHash=groupBySubstrate($featureArray);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;
  
  my $self={};
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#----------------------------------------------------------------
# $featureArray=$reader->loadGFF($filename);
#
sub loadGFF
  {
    my ($self,$gffFilename)=@_;
    my @features;
    if($gffFilename=~/\.gz$/)
      {open(GFF,"cat $gffFilename|gunzip|") || die $gffFilename}
    else
      {open(GFF,$gffFilename) || die $gffFilename}
    while(<GFF>)
      {
	next unless $_=~/\S+/;
	next if $_=~/^\s*\#/;
	my @fields=split/\s+/,$_;
	next unless @fields>7;
	my @additionalFields=splice(@fields,8);
	my $feature=new Feature($fields[0],$fields[1],$fields[2],
				$fields[3]-1,$fields[4],$fields[5],
				$fields[6],$fields[7],$fields[8],
				\@additionalFields);
	push @features,$feature;
      }
    close(GFF);
    @features=sort {$a->{begin} <=> $b->{begin}} @features;
    return \@features;
  }
#--------------------------------------------------------------------------
#   $bySubstrateHash=$reader->hashBySubstrate($filename);
sub hashBySubstrate
{
    my ($self,$filename)=@_;
    my $allFeatures=$self->loadGFF($filename);
    my $hash={};
    my $n=@$allFeatures;
    for(my $i=0 ; $i<$n ; ++$i)
    {
        my $feature=$allFeatures->[$i];
        my $substrate=$feature->getSubstrate();
        push @{$hash->{$substrate}},$feature;
    }
    return $hash;
}
#--------------------------------------------------------------------------
#   $groupArray=groupByOverlaps($featureArray); # each group=array of features
sub groupByOverlaps {
  my ($features)=@_;
  my $numFeatures=@$features;
  my $groups=[]; # each group is an array of GffFeature
  my $groupBounds=[]; # each element is an array [begin,end]
  my @features=sort {$a->getBegin() <=> $b->getBegin()} @$features;
  for(my $i=0 ; $i<$numFeatures ; ++$i) {
    my $feature=$features[$i];
    push @$groups,[$feature];
    push @$groupBounds,[$feature->getBegin(),$feature->getEnd()];
  }
  my $numGroups=@$groups;
  my $changes=1;
  while($changes) {
    $changes=0;
    for(my $i=0 ; $i+1<$numGroups ; ++$i) {
      my $thisGroup=$groups->[$i];
      my $nextGroup=$groups->[$i+1];
      if(groupsOverlap($groupBounds->[$i],$groupBounds->[$i+1])) {
	push @$thisGroup,@$nextGroup;
	my $bounds=$groupBounds->[$i];
	my $nextBounds=$groupBounds->[$i+1];
	$bounds->[0]=min($bounds->[0],$nextBounds->[0]);
	$bounds->[1]=max($bounds->[1],$nextBounds->[1]);
	splice(@$groups,$i+1,1);
	splice(@$groupBounds,$i+1,1);
	--$numGroups;
	$changes=1;
      }
    }

  }
  return $groups;
}
#--------------------------------------------------------------------------
sub groupsOverlap {
  my ($a,$b)=@_;
  return $a->[0]<$b->[1] && $b->[0]<$a->[1];
}
#--------------------------------------------------------------------------
sub min {
  my ($a,$b)=@_;
  return $a<$b ? $a : $b;
}
#--------------------------------------------------------------------------
sub max {
  my ($a,$b)=@_;
  return $a>$b ? $a : $b;
}
#--------------------------------------------------------------------------
#   $groupHash=groupBySubstrate($featureArray);
sub groupBySubstrate {
  my ($features)=@_;
  my $hash={};
  my $n=@$features;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $f=$features->[$i];
    push @{$hash->{$f->getSubstrate()}},$f;
  }
  return $hash;
}
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
#--------------------------------------------------------------------------

1;

