#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/genomics/perl','/home/bmajoros/perlib');
use Translation;
use FastaReader;
use FastaWriter;
use FileHandle;

my $usage="$0 <in.fasta> <out.fasta>";
die "$usage\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $fastaWriter=new FastaWriter;
my $filehandle=new FileHandle(">$outfile") || die "can't create $outfile";

my $fastaReader=new FastaReader($infile);
while(1)
{
    my ($defline,$seq)=$fastaReader->nextSequence();
    last unless defined $defline;
    my $orf;
    if($defline=~/startCodon=(\d+)/)
      {
	my $start=$1;
	$orf=substr($seq,$start,length($seq)-$start);
	#print "DEF=$defline\nRAW=$seq\nORF=$orf\n";
      }
    else {die;$orf=findLongestOrf(\$seq)}
    my $protein=Translation::translate(\$orf);
    if(defined $orf)
      {$fastaWriter->addToFasta($defline,$protein,$filehandle)}
    else
      {print STDERR "WARNING! No ORF found for $defline\n"}
}
close($filehandle);

#---------------------------------------------------------
sub findLongestOrf
{
    my ($seq)=@_;
    my $len=length $$seq;
    my ($longestOrf,$longestLength);
    for(my $frame=0 ; $frame<3 ; ++$frame)
    {
	my $orf=longestOrfInFrame(substr($$seq,$frame,$len-$frame));
	my $len=length $orf;
	if($len>$longestLength)
	{
	    $longestOrf=$orf;
	    $longestLength=$len;
	}
    }
    return $longestOrf;
}

sub longestOrfInFrame
{
    my ($seq)=@_;

    my ($open,$orfLen,$longestOrf,$longestLen,$orf);
    my $seqLen=length $seq;
    my $numCodons=int($seqLen/3);
    for(my $i=0 ; $i<$numCodons ; ++$i)
    {
	my $codon=substr($seq,3*$i,3);
	if($open)
	{
	    $orf.=$codon;
	    ++$orfLen;
	    if($codon eq "TAG" || $codon eq "TAA" || $codon eq "TGA")
	    {
		$open=0;
		if(!defined $longestLen || $orfLen>$longestLen)
		{
		    $longestLen=$orfLen;
		    $longestOrf=$orf;
		}
		undef $orf;
	    }
	}
	elsif($codon eq "ATG")
	{
	    $open=1;
	    $orfLen=1;
	    $orf.=$codon;
	}
    }
    return $longestOrf;
}
