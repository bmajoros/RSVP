#!/usr/bin/perl
use strict;
use SGE;

my $sge=new SGE();
my $n=$sge->countJobs("shufex");
print "$n shufex\n";
$n=$sge->countJobs("QRLOGIN");
print "$n QRLOGIN\n";
$n=$sge->countJobs("");
print "$n total\n";
