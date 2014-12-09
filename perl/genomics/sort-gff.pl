#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl');
use GffReader;

my $usage="$0 <*.gff>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my $reader=new GffReader;
my $array=$reader->loadGFF($filename);
my @sorted=sort 
  {
    $a->{fivePrime}!=$b->{fivePrime} ?
      $a->{fivePrime} <=> $b->{fivePrime} :
	$a->{threePrime} <=> $b->{threePrime}
  } @$array;

my $n=@sorted;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $feature=$sorted[$i];
    print $feature->toGff();
  }


