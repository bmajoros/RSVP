#!/usr/bin/perl
use strict;
use GffReader;

my $usage="$0 *.gff *.gff (subtracts the second from the first)";
die "$usage\n" unless @ARGV==2;
my ($firstFile,$secondFile)=@ARGV;

my ($vertices,$edges)=loadGff($secondFile);

if($firstFile=~/\.gz$/)
  {open(IN,"cat $firstFile|gunzip|") || die "can't open $firstFile"}
else
  {open(IN,$firstFile) || die "can't open $firstFile"}
while(<IN>)
  {
    next if(/^\s*\#/);
    my @line=split/\s+/,$_;
    my ($source,$featureType,$begin,$end,$strand)=
      ($line[1],$line[2],$line[3],$line[4],$line[6]);

    my $signature="$source($featureType:$begin-$end)";
    if($source=~/vertex/ && $vertices->{$signature}) 
      {print "===>FOUND: $signature\n";next}
    if($source=~/edge/ && $edges->{$signature}) 
      {print "===>FOUND: $signature\n";next}
    print "$signature\n";
  }
close(IN);

#-----------------------------------------------------------
sub loadGff
  {
    my ($filename)=@_;
    my (%vertices,%edges);

    if($filename=~/\.gz$/)
      {open(IN,"cat $filename|gunzip|") || die "can't open $filename"}
    else
      {open(IN,$filename) || die "can't open $filename"}
    while(<IN>)
      {
	next if(/^\s*\#/);
	my @line=split/\s+/,$_;
	my ($source,$featureType,$begin,$end,$strand)=
	  ($line[1],$line[2],$line[3],$line[4],$line[6]);

	### TEMPORARY HACK:
	if($featureType eq "initial-exon" || $featureType eq "single-exon")
	  {if($strand eq "+"){$begin-=3}else{$end+=3}}
	if($featureType eq "final-exon" || $featureType eq "single-exon")
	  {if($strand eq "+"){$end+=3}else{$begin-=3}}
	if($featureType eq "INTRON")
	  {$begin-=2;$end+=2}
	###

	my $signature="$source($featureType:$begin-$end)";
	if($source=~/vertex/) {$vertices{$signature}=1}
	elsif($source=~/edge/) {$edges{$signature}=1}
	else {print STDERR "source not recognized: $source\n"}
      }
    close(IN);
    return (\%vertices,\%edges);
  }
#-----------------------------------------------------------
sub old
  {
my $reader=new GffReader();
my $firstSet=$reader->loadGFF($firstFile);
my $secondSet=$reader->loadGFF($secondFile);

my (%vertices,%edges);
foreach my $feature (@$secondSet)
  {
    my $featureType=$feature->{featureType};
    if($featureType eq "edge")
      {
	my $source=$feature->{source};
	my $begin=$feature->{begin};
	my $end=$feature->{end};
	$edges{"$source:$begin-$end"}=1;
      }
    elsif($featureType eq "vertex")
      {
	my $source=$feature->{source};
	my $begin=$feature->{begin};
	my $end=$feature->{end};
	$vertices{"$source:$begin-$end"}=1;
      }
    elsif($featureType eq "transcript") {next}
    else
      {
	die "unknown feature type $featureType: ".$feature->toGff();
      }
  }

foreach my $feature (@$firstSet)
  {
    my $featureType=$feature->{featureType};
    if($featureType eq "edge")
      {
	my $source=$feature->{source};
	my $begin=$feature->{begin};
	my $end=$feature->{end};
	if(!defined($edges{"$source:$begin-$end"}))
	  {print $feature->toGff()}
      }
    elsif($featureType eq "vertex")
      {
	my $source=$feature->{source};
	my $begin=$feature->{begin};
	my $end=$feature->{end};
	if(!defined($vertices{"$source:$begin-$end"}))
	  {print $feature->toGff()}
      }
    elsif($featureType eq "transcript") {next}
  }
}


