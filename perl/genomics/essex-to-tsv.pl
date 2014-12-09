#!/usr/bin/perl
use strict;
use ProgramName;
use EssexParser;

my $name=ProgramName::get();
die "$name <in.sx> <field1> <field2> ...\n" unless @ARGV>=2;
my $infile=shift @ARGV;

my @fields=@ARGV;
my $numFields=@fields;
for(my $i=0 ; $i<$numFields ; ++$i) {
  my $field=$fields[$i];
  print "$field";
  if($i<$numFields-1) {print "\t"}
}
print "\n";

my $parser=new EssexParser($infile);
while(1) {
  my $root=$parser->nextElem();
  last unless $root;
  for(my $i=0 ; $i<$numFields ; ++$i) {
    my $tag=$fields[$i];
    my $elem=$root->findChild($tag);
    if(!$elem) {print "."}
    else {
      my $numChildren=$elem->numElements();
      for(my $j=0 ; $j<$numChildren ; ++$j) {
	my $val=$elem->getIthElem($j);
	if(EssexNode->isaNode($val)) {$val="*"}
	elsif($val eq "") {$val="."}
	print $val;
	if($j<$numChildren-1) {print ","}
      }
    }
    if($i<$numFields-1) {print "\t"}
  }
  print "\n";
}
