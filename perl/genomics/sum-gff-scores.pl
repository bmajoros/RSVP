#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl');
use GffReader;

my $usage="$0 <*.gff> <begin> <end>\n{both coordinates are 1-based and inclusive}\n";
die "$usage\n" unless @ARGV==3;
my ($filename,$begin,$end)=@ARGV;

my $sum;
my $reader=new GffReader;
my $array=$reader->loadGFF($filename);
my $n=@$array;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $feature=$array->[$i];
    if($feature->overlaps($begin,$end)) {$sum+=$feature->{score}}
  }
print "$sum\n";

