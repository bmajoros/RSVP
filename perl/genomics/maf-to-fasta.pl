#!/usr/bin/perl
use strict;
use ProgramName;
use FastaWriter;
use FileHandle;

my $name=ProgramName::get();
die "$name <in.maf>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $writer=new FastaWriter;
my $fh=new FileHandle($infile);
die "Can't open $infile\n" unless $fh;
while($_=<$fh>) {
  my @fields=split/\s+/,$_;
  if($fields[0] eq "s") {
    my $name=$fields[1];
    my $seq=$fields[6];
    my $defline=">$name";
    $writer->addToFasta($defline,$seq,\*STDOUT);
  }
  elsif($fields[0] eq "a") { print "\n" }
}


