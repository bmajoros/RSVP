#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <infile> <outfile>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

open(OUT,">$outfile") || die "can't write to $outfile\n";
open(IN,$infile) || die "can't open $infile\n";
<IN>;
print OUT "a score=1\n";
my $prevBlank=0;
while(<IN>) {
  chomp;
  if(/^\s*$/) {
    if(!$prevBlank) {print OUT "\n"}
    $prevBlank=1;
  }
  else { $prevBlank=0 }
  next if(/^[\s\*]*$/);
  if(/^\S/) {
    my @fields=split/\s+/,$_;
    die "can't parse line: \"$_\"" unless @fields==2;
    my $species=$fields[0];
    my $seq=$fields[1];
    my $len=length($seq);
    print OUT "s $species 0 $len + . $seq\n";
  }
}
close(IN);
close(OUT);



