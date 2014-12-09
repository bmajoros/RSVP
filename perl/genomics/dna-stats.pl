#!/usr/bin/perl
use strict;
use FastaReader;
use ProgramName;

my $name=ProgramName::get();
my $usage="$name <*.fasta>";
die "$usage\n" unless @ARGV==1;
my ($filename)=@ARGV;

my ($totalAs,$totalTs,$totalCs,$totalGs);
my $reader=new FastaReader($filename);
my $numSeqs;
my $lenSum;
while(1)
  {
    my ($defline,$sequence)=$reader->nextSequence();
    last unless defined $defline;

    my $numAs=($sequence=~s/A/A/g);
    my $numTs=($sequence=~s/T/T/g);
    my $numCs=($sequence=~s/C/C/g);
    my $numGs=($sequence=~s/G/G/g);
    my $total=$numAs+$numTs+$numCs+$numGs;
    next unless $total>0;
    ++$numSeqs;
    my $thisLen=length($sequence);
    $lenSum+=$thisLen;
    next unless $thisLen>0;

    $totalAs+=$numAs;
    $totalCs+=$numCs;
    $totalGs+=$numGs;
    $totalTs+=$numTs;

    my $percentA=int(100*$numAs/$total+5/9)/100;
    my $percentT=int(100*$numTs/$total+5/9)/100;
    my $percentC=int(100*$numCs/$total+5/9)/100;
    my $percentG=int(100*$numGs/$total+5/9)/100;
    
    my $pA=$numAs/$total;
    my $pT=$numTs/$total;
    my $pC=$numCs/$total;
    my $pG=$numGs/$total;
    
    my $H=entropy($pA,$pC,$pG,$pT);
    my $Hmax=-lg(1/4.0);
    my $Hnorm=$H/$Hmax;
    
    my $percentGC=$percentG+$percentC;
    my $percentAT=$percentA+$percentT;
    
    my $meanLen=int($lenSum/$numSeqs+5/9);
    
    $defline=~/^>(\S+)/;
    my $id=$1;
    print "$id L=$total A=$percentA T=$percentT C=$percentC G=$percentG GC=$percentGC AT=$percentAT H=$H H/Hmax=$Hnorm\n";
  }

#print "=====================================================\n";

my $total=$totalAs+$totalTs+$totalCs+$totalGs;
my $percentA=int(100*$totalAs/$total+5/9)/100;
my $percentT=int(100*$totalTs/$total+5/9)/100;
my $percentC=int(100*$totalCs/$total+5/9)/100;
my $percentG=int(100*$totalGs/$total+5/9)/100;

my $pA=$totalAs/$total;
my $pT=$totalTs/$total;
my $pC=$totalCs/$total;
my $pG=$totalGs/$total;

my $H=entropy($pA,$pC,$pG,$pT);
my $Hmax=-lg(1/4.0);
my $Hnorm=$H/$Hmax;

my $percentGC=$percentG+$percentC;
my $percentAT=$percentA+$percentT;

my $meanLen=int($lenSum/$numSeqs+5/9);

print "AVERAGE L=$meanLen A=$percentA T=$percentT C=$percentC G=$percentG GC=$percentGC AT=$percentAT H=$H H/Hmax=$Hnorm N=$numSeqs\n";


sub entropy
  {
    my @array=@_;
    my $H=0.0;
    foreach my $p (@array) {$H+=entropyTerm($p)}
    return $H;
  }

sub entropyTerm
  {
    my ($p)=@_;
    return $p>0 ? -$p*lg($p) : 0.0;
  }

sub lg
  {
    my ($x)=@_;
    return log($x)/log(2);
  }
