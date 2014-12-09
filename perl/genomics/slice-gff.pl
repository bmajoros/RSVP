#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use GffReader;

my $usage="$0 <*.gff> <begin> <end>\n{begin & end are one/base-based}";
die "$usage\n" unless @ARGV==3;
my ($filename,$begin,$end)=@ARGV;

my $reader=new GffReader;
my $features=$reader->loadGFF($filename);
my $n=@$features;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $feature=$features->[$i];
    my $trimmed=$feature->intersect($begin,$end);
    next unless defined $trimmed;
    print $trimmed->toGff();
  }

