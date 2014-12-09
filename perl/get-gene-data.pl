#!/usr/bin/perl
#Parses a list of genes (in GFF format), selects the first n single-isoform
#genes not on Chr1, ChrM, or ChrC, and creates pileup and GFF files for those
#genes. Also creates a GFF file with a list of all the selected genes.

#The output gene list file is currently called single-transcript-genes.gff,
#but the name can be changed by modifying the $outGeneListName variable.
#The skipped chromosomes can also be changed by modifying the if statement
#labeled with "HARD-CODED SKIP CHROMOSOMES."
#
#Currently the script only generates pileup and GFF files for the genes. To
#output SAM files as well, uncomment the system() call at the bottom of the
#script. Note that these SAM files will contain reads from both strands. 
#
#As of 6/6/12, this script selects only high-coverage genes
#for the histogram training set, using code copied from check-coverage.pl.

use warnings;
use strict;
use GffTranscriptReader;
use Getopt::Std;

our $opt_s;
getopts('s:');

die "usage: get-gene-data.pl [-s samtoolspath] <output-folder> <input-gff-file> <bam-file> num-genes threshold" unless @ARGV == 5;

my $samtoolsPath = "samtools";
if (defined($opt_s))
{
  $samtoolsPath = "$opt_s/samtools";
}


my $outDir = $ARGV[0];
my $inGeneList = $ARGV[1];
my $bamFile = $ARGV[2];
my $numGenes = $ARGV[3];
my $threshold = $ARGV[4];
my $outGeneListName = "single-transcript-genes.gff";

my $reader = new GffTranscriptReader();
$reader -> setStopCodons({TAG=>1,TAA=>1,TGA=>1});

my %chrHash = %{$reader -> hashBySubstrate("$inGeneList")};
my @chromosomes = keys(%chrHash);

my @singleTranscripts;

for my $key (@chromosomes)
{
    if ($key eq "ChrM" || $key eq "Chr1" || $key eq "ChrC") #HARD-CODED SKIP CHROMOSOMES
    {
	next;
    }
    
    my @trRefArray = @{$chrHash{$key}};

    for my $trRef (@trRefArray)
    {
	my $geneRef = $trRef -> getGene();
	my $numTranscripts = $geneRef -> getNumTranscripts();
    
	if ($numTranscripts == 1)
	{
	    push(@singleTranscripts, $geneRef);
	}
	my $id = $geneRef -> getId();    
    }

}

open(GENE_GFF_FILE, "> $outDir/$outGeneListName") || die ("Can't open gene output file");

my $chunkId = 1;

for (my $counter = 0; $counter < $numGenes; $counter++)
{
    if ($counter >= @singleTranscripts)
    {
        last;
    }

    my $sGeneRef = $singleTranscripts[$counter];
    if ($chunkId > $numGenes)
    {
	last;
    }
    my $transcript = $sGeneRef -> getIthTranscript(0);
    my $transcriptID = $transcript -> getTranscriptId();
    my $begin = $transcript -> getBegin();
    my $end = $transcript -> getEnd();
    my $strand = $transcript -> getStrand();
    (my $chrNum) = ($transcriptID =~ /AT(\w)G/);
    my $chromosome = "Chr$chrNum";
    my $outputString = "Chr$chrNum\tarab\ttraining\t$begin\t$end\t.\t$strand\t.\ttranscript_id=$transcriptID\tchunk_id=$chunkId;\n";
    
    my $gffString = $transcript -> toGff();

    #print statements at bottom of loop

    #check for overlapping genes

    my $nGeneRef = $singleTranscripts[$counter + 1];
    my $nTranscript;
    my $nBegin;
    my $nStrand;
    if (defined $nGeneRef)
    {
	$nTranscript = $nGeneRef -> getIthTranscript(0);
	$nBegin = $nTranscript -> getBegin();
	$nStrand = $nTranscript -> getStrand();
    }

    my $pGeneRef = $singleTranscripts[$counter - 1];
    my $pTranscript;
    my $pEnd;
    my $pStrand;
    if (defined $pGeneRef)
    {
	$pTranscript = $pGeneRef -> getIthTranscript(0);
	$pEnd = $pTranscript -> getEnd();
	$pStrand = $pTranscript -> getStrand();
    }

    if (defined $nBegin)
    {
	if ($nBegin < $end)
	{
	    $numGenes++; #so we still have the proper number of genes
	    next;
	}
    }
    if (defined $pEnd)
    {
	if ($pEnd > $begin)
	{
	    $numGenes++;
	    next;
	}
    }

    #check coverage
    my @sam = `$samtoolsPath view -h $bamFile '$chromosome:$begin-$end'`;
    my $numReads = 0;
    for (my $i = 0; $i < @sam; $i++)
    {
        my $line = $sam[$i];
        my @lineSplit = split(/\s+/, $line);
        my $coord = $lineSplit[3];
        if ($line =~ /^@/) {next;}
        $numReads++;
    }
    my $length = $end - $begin + 1;
    my $rpkb = $numReads / ($length / 1000);
    if($rpkb < $threshold)
    {
        $numGenes++;
        next;
    }
 
    #adjust 3' buffer size based on neighboring genes
    
    my $bufferBegin = $begin;
    my $bufferEnd = $end;
    
    if ($strand eq "+")
    {
	$bufferEnd = $end + 4000;
	if (defined $nBegin)
	{
	    my $tempEnd;
	    
	    $tempEnd = int(($nBegin - $end) / 2) + $end; #truncation instead of rounding to make sure buffers don't overlap

	    if ($tempEnd < $bufferEnd)
	    {
		$bufferEnd = $tempEnd;
	    }
	}
	if ($bufferEnd - $end > 4000) #debug
	{
	    print ("flag! chunk $chunkId\n");
	}
    }
    else
    {
	$bufferBegin = $begin - 4000;
	if (defined $pEnd)
	{
	    my $tempBegin;
	    
	    $tempBegin = int(($begin - $pEnd) / 2) + 1 + $pEnd;
	    
	    if ($tempBegin > $bufferBegin)
	    {
		$bufferBegin = $tempBegin;
	    }
	}
    }

    print (GENE_GFF_FILE $outputString);

    open(EXON_GFF_FILE, "> $outDir/$chunkId.gff");
    print (EXON_GFF_FILE $gffString);

    system("$samtoolsPath mpileup -r '$chromosome:$bufferBegin-$bufferEnd' $bamFile > $outDir/$chunkId.pileup");
	
    #system("$samtoolsPath view -h $bamFile $chromosome:$bufferBegin-$bufferEnd > $outDir/$chunkId.sam");
    $chunkId++;
}
