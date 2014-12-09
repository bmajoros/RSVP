#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use ProgramName;

use Getopt::Std;
our ($opt_i,$opt_d,$opt_s);
getopts('ids');

my $name=ProgramName::get();
my $usage="

$name [-ids] \"pattern\" <*.fasta> D/S
  D = apply pattern to defline
  S = apply pattern to sequence
  -i = report sequence index only
  -d = report deflines only
  -s = report sequences only
";
die "$usage\n" unless @ARGV==3;
my ($pattern,$filename,$code)=@ARGV;
die "$usage\n" unless($code eq "D" || $code eq "S");

my $writer=new FastaWriter();
my $reader=new FastaReader($filename);
$reader->dontUppercase();
my $seqNum=-1;
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;
    ++$seqNum;
    if($code eq "D")
      {next unless $defline=~/$pattern/}
    else
      {next unless $sequence=~/$pattern/}
    if($opt_d) {print "$defline"}
    elsif($opt_s) {print "$sequence\n"}
    elsif($opt_i) {print "$seqNum\n"}
    else {
      $writer->addToFasta($defline,$sequence,\*STDOUT);
    }
  }


