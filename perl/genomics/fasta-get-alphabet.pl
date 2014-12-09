#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use FastaReader;

my $usage="$0 <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my (%hash,$total);
open(IN,$filename) || die;
while(<IN>)
  {
    chop;
    unless(/>/)
      {
	my @a=split//,$_;
	foreach my $x (@a) {++$hash{$x};++$total}
      }
  }
close(IN);
my @a=keys %hash;
my $alphabet=join('',@a);
print "$alphabet\n\n";
foreach my $x (@a)
  {
    my $count=$hash{$x};
    my $p=int(1000*$count/$total)/10;
    print "$x $p%\t($count)\n";
  }

