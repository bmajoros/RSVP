#!/usr/bin/perl
use strict;
use ProgramName;
use Getopt::Std;
our $opt_b;
getopts('b');

my $name=ProgramName::get();
die "
cat datafile | $name [-b] mapfile > outfile
    -b = keep BOTH the old and new identifiers

" unless @ARGV==1;
my ($mapfile)=@ARGV;

# LOAD THE MAP ===========================================================
my %map;
open(IN,$mapfile) || die "can't open $mapfile";
while(<IN>)
  {
    my @fields=split/\s+/,$_;
    next unless @fields==2;
    $map{$fields[0]}=$fields[1];
  }
close(IN);

# PROCESS THE INPUT STREAM ===============================================
while(<STDIN>)
  {
    if(/^\s*(\S+)(.*)/)
      {
	my ($id,$rest)=($1,$2);
	my $newID=$map{$id};
	next unless $newID;
	if($opt_b) 
	  {print "$newID\t$id$rest\n"}
	else
	  {
	    $rest=~s/^\s+//;
	    print "$newID\t$rest\n"
	  }
      }
  }









