#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use TempFilename;

my $usage="$0 <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $tempFile=TempFilename::generate();
open(OUT,">$tempFile") || die "Can't create temporary file \"$tempFile\"\n";
my $writer=new FastaWriter;

my $hash=FastaReader::readAllAndKeepDefs($infile);
my @keys=keys %$hash;
my $n=@keys;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $key=$keys[$i];
    my $entry=$hash->{$key};
    my ($def,$seq)=@$entry;
    $writer->addToFasta($def,$seq,\*OUT);
  }

close(OUT);

system("mv $tempFile $infile");

