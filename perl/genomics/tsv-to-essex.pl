#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.tsv> <top-level-tag>\n" unless @ARGV==2;
my ($infile,$rootTag)=@ARGV;

open(IN,$infile) || die "can't open file \"$infile\"\n";
my $headerLine=<IN>;
chomp $headerLine;
my @fields=split/\t/,$headerLine;
my $numFields=@fields;
while(<IN>) {
  chomp;
  my @data=split/\t/,$_;
  print "($rootTag\n";
  for(my $i=0 ; $i<$numFields ; ++$i) {
    my $field=$fields[$i];
    if($i==$numFields-1 && @data>$numFields) {
      print "    ($field";
      for(my $j=$i ; $j<@data ; ++$j) {
	my $datum=$data[$j];
	print " $datum";
      }
      print ")\n";
    }
    else {
      my $datum=$data[$i];
      print "    ($field $datum)\n";
    }	
  }
  print ")\n";
}
close(IN);


