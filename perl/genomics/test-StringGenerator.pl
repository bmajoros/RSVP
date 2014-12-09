#!/usr/bin/perl -w
use strict;
use StringGenerator;

#my $g=new StringGenerator({numSymbols=>6});
my $g=new StringGenerator({distribution => 
			   {
			    a=>.1,
			    b=>.2,
			    c=>.3,
			    d=>.15,
			    e=>.05,
			    f=>.2
			   }});

my $s=$g->generateString(10000);

my $n;
my $p;
my $len=length $$s;
$n=($$s=~s/a/a/g); $p=$n/$len; print "a $p\n";
$n=($$s=~s/b/b/g); $p=$n/$len; print "b $p\n";
$n=($$s=~s/c/c/g); $p=$n/$len; print "c $p\n";
$n=($$s=~s/d/d/g); $p=$n/$len; print "d $p\n";
$n=($$s=~s/e/e/g); $p=$n/$len; print "e $p\n";
$n=($$s=~s/f/f/g); $p=$n/$len; print "f $p\n";

