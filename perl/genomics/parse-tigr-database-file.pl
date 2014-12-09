#!/usr/bin/perl
use strict;

$0=~/([^\/]+)\s*$/;
my $usage="$1 <infile>";
die "$usage\n" unless @ARGV==1;
my ($infile)=@ARGV;

open(IN,$infile) || die "can't open $infile\n";
my @fields;
while(<IN>)
  {
    if(/^[\s-]*-[\s-]*$/) {last}
    if(/^\s*$/) {push @fields,["",0,0]}
    my $line=$_;
    while($line=~/^(\s*)(\S+\s*)(.*)/)
      {
	my ($spacer,$field,$remainder)=($1,$2,$3);
	my $fieldLength=length $field;
	$line=$remainder;
	$field=~s/\s+//g;
	my $fieldBegin=length($spacer);
	print "FIELD=$field LENGTH=$fieldLength BEGIN=$fieldBegin SPACER=\"$spacer\"\n";
	push @fields,[$field,$fieldLength,$fieldBegin];
      }
  }
my $recordNum=0;
while(<IN>)
  {
    if(/^[\s-]+$/) {next}
    print "\nrecord #$recordNum\n";
    my $recordText;
    foreach my $field (@fields)
      {
	if(length($_)<2) {$_=<IN>}
	my ($label,$length,$begin)=@$field;
	if($label eq "") 
	  {
	    $_=~s/\s*$//g;
	    chomp $recordText;
	    $recordText.="$_\n";
	    $_=<IN>; 
	    next
	  }
	my $data=substr($_,$begin,$length);
	$_=substr($_,$begin+$length,length($_)-$length-$begin);
	chomp $data;
	$data=~s/^\s*//g;
	#print "$label=$data\n";
	$recordText.="$label=$data\n";
      }
    print "$recordText";
    ++$recordNum;
  }
close(IN);




