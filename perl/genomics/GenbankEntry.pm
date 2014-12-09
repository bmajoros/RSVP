package GenbankEntry;
use strict;

######################################################################
#
# GenbankEntry.pm bmajorostigr.org 11/16/2004
#
# Just an array of (key,value) pairs, where each value is either a
# string, a number, or a nested GenbankEntry.  NOTE that the keys are
# *not* assumed to be unique!
#
#
# Attributes:
#
# Methods:
#   $entry=new GenbankEntry();
#   $entry->addPair($key,$value);
#   my $n=$entry->numPairs();
#   my $pair=$entry->getIthPair($i); # returns pointer to [key,value]
#   my $pairs=$entry->findPairs($key); # returns array of pairs
#   my $value=$entry->findUnique($key); # returns one value element
#   $entry->print(*STDOUT);
#   my $cds=$entry->getCDS(); # returns array of [begin,end] pairs
#   my $cds=$entry->getExons(); # returns array of [begin,end] pairs
#   my $seq=$entry->getSubstrate();
#   $entry->consolidateFeatures(); # combines FEATURES & COMMENTS elements
#
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#   $entry=new GenbankEntry();
sub new
{
  my ($class)=@_;

  my $self=
    {
     array=>[]
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $entry->addPair($key,$value);
sub addPair
  {
    my ($self,$key,$value)=@_;
    my $array=$self->{array};
    push @$array,[$key,$value];
  }
#---------------------------------------------------------------------
#   my $n=$entry->numPairs();
sub numPairs
  {
    my ($self)=@_;
    my $array=$self->{array};
    my $n=@$array;
    return $n;
  }
#---------------------------------------------------------------------
#   my $pair=$entry->getIthPair($i); # returns pointer to [key,value]
sub getIthPair
  {
    my ($self,$i)=@_;
    my $array=$self->{array};
    return $array->[$i];
  }
#---------------------------------------------------------------------
#   my $pairs=$entry->findPairs($key); # returns array of pairs
sub findPairs
  {
    my ($self,$key)=@_;
    my $array=$self->{array};
    my $n=@$array;
    my $pairs=[];
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	if($pair->[0] eq $key) {push @$pairs,$pair}
      }
    return $pairs;
  }
#---------------------------------------------------------------------
#   my $value=$entry->findUnique($key); # returns one value element
sub findUnique
  {
    my ($self,$key)=@_;
    my $array=$self->{array};
    my $n=@$array;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	if($pair->[0] eq $key) {return $pair->[1]}
      }
    return undef;
  }
#---------------------------------------------------------------------
sub print
  {
    my ($self,$handle,$level)=@_;
    my $array=$self->{array};
    my $n=@$array;
    my $pad=' 'x($level*3);
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $pair=$array->[$i];
	my ($key,$value)=@$pair;
	print "$pad$key => ";
	if(UNIVERSAL::isa($value,"GenbankEntry"))
	  {
	    print "\n$pad   \{\n";
	    $value->print($handle,$level+1);
	    print "$pad   }\n";
	  }
	else {print "$value\n"}
      }
  }
#---------------------------------------------------------------------
#   my $cds=$entry->getCDS(); # returns array of [begin,end] pairs
sub getCDS
  {
    my ($self)=@_;
    my $features=$self->findUnique("FEATURES");
    #die "No FEATURES element found in Genbank entry!\n" unless $features;
    if(!$features) {return undef}
    my $seq=$self->getSubstrate();
    my $len=length($seq);
    my $rhs=$features->findUnique("CDS");
    return undef unless $rhs;
    #die "No CDS element found in Genbank FEATURES clause!\n" unless $rhs;
    my $array=[];
    if($rhs=~/^\s*(complement\()?([A-Za-z0-9]+:)?(<?)(\d+)\.\.(>?)(\d+)/)
      {
	my ($complement,$junk,$less,$begin,$greater,$end)=
	  ($1,$2,$3,$4,$5,$6);
	if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	  {
	    push @$array,[1,$len];
	    return $array;
	  }
	if($complement) {($begin,$end)=($end,$begin)}
	push @$array,[$begin,$end];
      }
    elsif($rhs=~/^\s*(complement\()?join\(([A-Za-z0-9]+:)?([^)]+)\)/)
      {
	my ($complement,$junk,$coords)=($1,$2,$3);
	if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	  {
	    push @$array,[1,$len];
	    return $array;
	  }
	$coords=~s/\s+//g;
	my @exons=split/,/,$coords;
	foreach my $exon (@exons)
	  {
	    if($exon=~/(<?)(\d+)\.\.(>?)(\d+)/)
	      {
		my ($less,$begin,$greater,$end)=($1,$2,$3,$4);
		if($complement) {($begin,$end)=($end,$begin)}
		push @$array,[$begin,$end];
	      }
	      elsif($exon=~/^\s*(\d+)\s*$/)
		{
		  my ($begin,$end)=($1,$1);
		  if($complement) {($begin,$end)=($end,$begin)}
		  push @$array,[$begin,$end];
		}
	      else {next}#die "Can't parse Genbank exon: $exon\n"}
	  }
      }
    else
      {die "Can't parse CDS clause in Genbank entry: $rhs\n"}
    return $array;
  }
#---------------------------------------------------------------------
#   my $cds=$entry->getExons(); # returns array of [begin,end] pairs
sub getExons
  {
    my ($self)=@_;
    my $features=$self->findUnique("FEATURES");
    if(!$features) {return undef}
    my $seq=$self->getSubstrate();
    my $len=length($seq);
    my $pairs=$features->findPairs("exon");
    my $morePairs=$features->findPairs("mRNA");
    push @$pairs,@$morePairs;
    my $array=[];
    foreach my $pair (@$pairs)
      {
	my $rhs=$pair->[1];
	return undef unless $rhs;
	if($rhs=~
	   /^\s*(complement\()?([A-Za-z0-9]+:)?(<?)(\d+)\.\.(>?)(\d+)/)
	  {
	    my ($complement,$junk,$less,$begin,$greater,$end)=
	      ($1,$2,$3,$4,$5,$6);
	    if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	      {
		push @$array,[1,$len];
		return $array;
	      }
	    if($complement) {($begin,$end)=($end,$begin)}
	    push @$array,[$begin,$end];
	  }
	elsif($rhs=~/^\s*(complement\()?join\(([A-Za-z0-9]+:)?([^)]+)\)/)
	  {
	    my ($complement,$junk,$coords)=($1,$2,$3);
	    if($junk)### TOTALLY BOGUS...FOR BURGE'S EST FILE
	      {
		push @$array,[1,$len];
		return $array;
	      }
	    $coords=~s/\s+//g;
	    my @exons=split/,/,$coords;
	    foreach my $exon (@exons)
	      {
		if($exon=~/(<?)(\d+)\.\.(>?)(\d+)/)
		  {
		    my ($less,$begin,$greater,$end)=($1,$2,$3,$4);
		    if($complement) {($begin,$end)=($end,$begin)}
		    push @$array,[$begin,$end];
		  }
		elsif($exon=~/^\s*(\d+)\s*$/)
		  {
		    my ($begin,$end)=($1,$1);
		    if($complement) {($begin,$end)=($end,$begin)}
		    push @$array,[$begin,$end];
		  }
		else {next}#die "Can't parse Genbank exon: $exon\n"}
	      }
	  }
	else {return undef}
	#{die "Can't parse exon clause in Genbank entry: $rhs\n"}
      }
    return $array;
  }
#---------------------------------------------------------------------
#   my $seq=$entry->getSubstrate();
sub getSubstrate
  {
    my ($self)=@_;
    my $seq=$self->findUnique("ORIGIN");
    $seq="\U$seq";
    return $seq;
  }
#---------------------------------------------------------------------
#   $entry->consolidateFeatures(); # combines FEATURES & COMMENTS elements
sub consolidateFeatures
  {
    my ($self)=@_;
    my $array=$self->{array};
    my $n=@$array;
    my $features;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $elem=$array->[$i];
	my ($key,$value)=@$elem;
	if(($key eq "FEATURES") || ($key eq "COMMENT"))
	  {
	    #print "$key\n";
	    if(defined $features)
	      {
		#print "$features\n";
		$features->subsume($value);
		$elem->[0]="FEATURES";
		splice(@$array,$i,0);
		--$n;
		--$i;
	      }
	    else {$features=$value}
	  }
      }
  }
#---------------------------------------------------------------------
sub subsume
  {
    my ($self,$other)=@_;
    my $array=$self->{array};
    my $otherArray=$other->{array};
    push @$array,@$otherArray;
  }
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------

1;

