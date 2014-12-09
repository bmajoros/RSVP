#!/usr/bin/perl
use strict;
use SummaryStats;
use FileSize;
use FastaReader;
use Numbers;
use GffTranscriptReader;

###############################################################
# This program extracts the following statistics from an entire
# genome: total coding length, total intron length, total
# intergenic length, total noncoding length, number of exons
# per gene, mean gene length, mean exon length, mean desert
# size
###############################################################

my $usage="$0 <gff-directory> <fasta-directory>";
die "$usage\n" unless @ARGV==2;
my ($gffDirectory,$fastaDirectory)=@ARGV;

my $gffFiles=`ls $gffDirectory/*.gff`;
my @gffFiles=split/\s+/,$gffFiles;
if($fastaDirectory=~/\/\s*$/) {chop $fastaDirectory}

my ($totalCoding,$totalIntron,$totalDesert,$totalNoncoding,
    $totalGeneExtent,$totalGenomeSize,@exonsPerGene);
my (@exonSizes,@intronSizes,@desertSizes,@geneSizes);
my $gffReader=new GffTranscriptReader;
my $totalTranscripts=0;
foreach my $gffFile (@gffFiles)
  {
    print STDERR "Loading transcript coordinates from $gffFile...\n";
    my $transcripts=$gffReader->loadGFF($gffFile);

    $gffFile=~/([^\/]+)\.gff/;
    my $filestem=$1;
    my $fastaFile="$fastaDirectory/$filestem.fasta";
    print STDERR "Processing $fastaFile...\n";
    die unless -e $fastaFile;
    my $axisLength=FastaReader::getGenomicSize($fastaFile);
    $totalGenomeSize+=$axisLength;

    print "$filestem axis length=$axisLength\n";
    my $numTranscripts=@$transcripts;
    $totalTranscripts+=$numTranscripts;
  TRANSCRIPT:
    for(my $j=0 ; $j<$numTranscripts ; ++$j)
      {
	my $transcript=$transcripts->[$j];
	if($j<$numTranscripts-1)
	  {
	    my $nextTranscript=$transcripts->[$j+1];

	    if($transcript->overlaps($nextTranscript)) {next TRANSCRIPT}###

	    if(!$transcript->overlaps($nextTranscript))
	      {
		my $desertSize=$nextTranscript->{fivePrime}-
		  $transcript->{threePrime};
		if($desertSize<0) {die "$transcript->{transcriptId}:$transcript->{fivePrime}-$transcript->{threePrime} is followed by $nextTranscript->{transcriptId}:$nextTranscript->{fivePrime}-$nextTranscript->{threePrime}"}
		push @desertSizes,$desertSize;
	      }
	  }
	my $strand=$transcript->{strand};
	my $exons=$transcript->{exons};
	my $numExons=@$exons;
	my $transcriptExtent=
	  $transcript->{threePrime}-$transcript->{fivePrime};
	my $transcriptLength=0;
	my $numExons=@$exons;
	push @exonsPerGene,$numExons;
	my $geneLen;
	for(my $i=0 ; $i<$numExons ; ++$i)
	  {
	    my $exon=$exons->[$i];
	    if($i<$numExons-1)
	      {
		my $nextExon=$exons->[$i+1];
		next if $nextExon->overlaps($exon);
		my $intronSize=
		  ($strand eq "+" ?
		  $nextExon->{begin}-$exon->{end} :
		  $exon->{begin}-$nextExon->{end});
		die "$transcript->{transcriptId} has an intron=$intronSize ($exon->{begin} minus $nextExon->{end}" unless $intronSize>=0;
		push @intronSizes,$intronSize;
	      }
	    my $exonLen=$exon->{end}-$exon->{begin};
	    $geneLen+=$exonLen;
	    die unless $exonLen>0;
	    $totalCoding+=$exonLen;
	    $transcriptLength+=$exonLen;
	    push @exonSizes,$exonLen;
	  }
	push @geneSizes,$geneLen;
	my $intronSum=$transcriptExtent-$transcriptLength;
	$totalIntron+=$intronSum;
	$totalGeneExtent+=$transcriptExtent;
	$totalNoncoding+=$intronSum;
      }
  }
$totalDesert=$totalGenomeSize-$totalGeneExtent;
$totalNoncoding=$totalIntron+$totalDesert;

my $desertPercent=int(100*$totalDesert/$totalGenomeSize+5/9);
my $codingPercent=int(100*$totalCoding/$totalGenomeSize+5/9);
my $noncodingPercent=int(100*$totalNoncoding/$totalGenomeSize+5/9);
my $extentPercent=int(100*$totalGeneExtent/$totalGenomeSize+5/9);
my $intronPercent=int(100*$totalIntron/$totalGenomeSize+5/9);

$totalNoncoding=Numbers::addCommas($totalNoncoding);
$totalDesert=Numbers::addCommas($totalDesert);
$totalGenomeSize=Numbers::addCommas($totalGenomeSize);
$totalCoding=Numbers::addCommas($totalCoding);
$totalIntron=Numbers::addCommas($totalIntron);
$totalGeneExtent=Numbers::addCommas($totalGeneExtent);
print "total genome length:     $totalGenomeSize bp\n";
print "total coding length:     $totalCoding bp\t($codingPercent\%)\n";
print "total gene extent:       $totalGeneExtent bp\t($extentPercent\%)\n";
print "total intron length:     $totalIntron bp\t($intronPercent\%)\n";
print "total intergenic length: $totalDesert bp\t($desertPercent\%)\n";
print "total noncoding length:  $totalNoncoding bp\t($noncodingPercent\%)\n";

my ($meanExon,$sdExon,$minExon,$maxExon)=
  SummaryStats::summaryStats(\@exonSizes);
$meanExon=int(100*$meanExon+5/9)/100;
$sdExon=int(100*$sdExon+5/9)/100;
my $n=@exonSizes;
print "Exon sizes: $meanExon+/-$sdExon ($minExon-$maxExon) N=$n\n";

my ($meanIntron,$sdIntron,$minIntron,$maxIntron)=
  SummaryStats::summaryStats(\@intronSizes);
$meanIntron=int(100*$meanIntron+5/9)/100;
$sdIntron=int(100*$sdIntron+5/9)/100;
my $n=@intronSizes;
print "Intron sizes: $meanIntron+/-$sdIntron ($minIntron-$maxIntron) N=$n\n";

my ($meanDesert,$sdDesert,$minDesert,$maxDesert)=
  SummaryStats::summaryStats(\@desertSizes);
$meanDesert=int(100*$meanDesert+5/9)/100;
$sdDesert=int(100*$sdDesert+5/9)/100;
my $n=@desertSizes;
print "Desert sizes: $meanDesert+/-$sdDesert ($minDesert-$maxDesert) N=$n\n";

my ($meanExonCount,$sdExonCount,$minExonCount,$maxExonCount)=
  SummaryStats::summaryStats(\@exonsPerGene);
$meanExonCount=int(100*$meanExonCount+5/9)/100;
$sdExonCount=int(100*$sdExonCount+5/9)/100;
my $n=@exonsPerGene;
print "#exons/gene: $meanExonCount+/-$sdExonCount ($minExonCount-$maxExonCount) N=$n\n";

my ($meanGeneSize,$sdGeneSize,$minGeneSize,$maxGeneSize)=
  SummaryStats::summaryStats(\@geneSizes);
$meanGeneSize=int(100*$meanGeneSize+5/9)/100;
$sdGeneSize=int(100*$sdGeneSize+5/9)/100;
my $n=@geneSizes;
print "Gene length (sum of intron lengths): $meanGeneSize+/-$sdGeneSize ($minGeneSize-$maxGeneSize) N=$n\n";

print "Total #transcripts: $totalTranscripts\n";

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------








