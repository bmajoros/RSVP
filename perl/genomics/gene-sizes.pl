#!/usr/bin/perl
use strict;
#use lib('/home/bmajoros/perlib','/home/bmajoros/genomics/perl');
use SummaryStats;
use GffTranscriptReader;

my $usage="$0 <*.gff> <*.fasta>";
die "$usage\n" unless @ARGV==2;
my ($gffFilename,$seqFilename)=@ARGV;

# Read the feature coordinates from the GFF file
my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($gffFilename);

my $numTranscripts=@$transcripts;
my (@numExons);
for(my $i=0 ; $i<$numTranscripts ; ++$i)
{
    my $transcript=$transcripts->[$i];
    my $exons=$transcript->{exons};
    my $transcriptId=$transcript->{transcriptId};
    my $numExons=@$exons;
    push @numExons,$numExons;
    my $threePrime=$transcript->{threePrime};
    my $fivePrime=$transcript->{fivePrime};
    my $extent=$threePrime-$fivePrime;
    my $strand=$transcript->{strand};
	
    my (@exonLengths,@intronLengths,$exonOutput,$transcriptLength);
    for(my $i=0 ; $i<$numExons ; ++$i)
    {
	my $exon=$exons->[$i];
	my $begin=$exon->{begin};
	my $end=$exon->{end};
	my $exonLen=$end-$begin;
	$transcriptLength+=$exonLen;
	push @exonLengths,$exonLen;
	$exonOutput.="\t$begin-$end $exonLen ";
	if($i<$numExons-1)
	{
	    my $nextExon=$exons->[$i+1];
	    my $intronSize=$nextExon->{begin}-$end;
	    $exonOutput.="$intronSize\n";
	    push @intronLengths,$intronSize;
	}
    }
    print "transcript \#$transcriptId ($strand) $numExons exons, extent=$extent len=$transcriptLength\n";

    my ($meanExonLen,$sdExonLen,$minExonLen,$maxExonLen)=
      SummaryStats::summaryStats(\@exonLengths);
    $meanExonLen=int($meanExonLen);
    $sdExonLen=int($sdExonLen);
    print "\texons: $meanExonLen+/-$sdExonLen ($minExonLen-$maxExonLen) ";

    if($numExons>1)
    {
	my ($meanIntronLen,$sdIntronLen,$minIntronLen,$maxIntronLen)=
	  SummaryStats::summaryStats(\@intronLengths);
	$meanIntronLen=int($meanIntronLen);
	$sdIntronLen=int($sdIntronLen);
	print "introns: $meanIntronLen+/-$sdIntronLen ($minIntronLen-$maxIntronLen)";
    }
    print "\n$exonOutput\n";
}



