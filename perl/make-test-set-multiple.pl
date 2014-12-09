#!/usr/bin/perl
#Modified version of Bill's script that use genes, not transcripts,
#to select the test set, and chooses only multiple-isoform genes.
use strict;
use FastaReader;
use CompiledFasta;
use FastaWriter;
use GffTranscriptReader;
use Translation;
use FileSize;
use ProgramName;
use Getopt::Std;
$|=1;

our $opt_c;
getopts('c:');

my $name=ProgramName::get();
my $usage="$name [-c substrate_id] <*.gff> <*.fasta> <max#> <margin-size> TAG,TGA,TAA

where -c indicates that the fasta file is compiled (*.cof)
         (compiled fasta files contain no defline and no whitespace)
";
die "$usage\n" unless @ARGV==5;
my ($gffFile,$fastaFile,$maxN,$margin,$stops)=@ARGV;
die "$gffFile does not exist\n" unless -e $gffFile;
die "$fastaFile does not exist\n" unless -e $fastaFile;

if(!-e "chunks") {system("mkdir chunks")}

my @stops=split/\s+/,$stops;
my %stops;
foreach my $stop (@stops) {$stops{$stop}=1}

##################### LOAD GFF #######################
my $gffReader=new GffTranscriptReader;
$gffReader->setStopCodons(\%stops);
my $genes=$gffReader->loadGenes($gffFile);
my $numGenes=@$genes;

# Index by contig ID and add padding around each gene
my %records;
for(my $i=0 ; $i<$numGenes ; ++$i)
  {
    my $gene=$genes->[$i];
    my $begin=$gene->getBegin();
    my $end=$gene->getEnd();
    die unless $begin<$end;
    next unless $begin>1;
    $begin-=$margin;
    $end+=$margin;
    my $substrate=$gene->getSubstrate();
    if(defined($opt_c)) {next unless $substrate eq $opt_c}
    my $strand=$gene->getStrand();
    push @{$records{"\L$substrate"}},[$begin,$end-$begin,$strand,$gene];
  }

# Eliminate any gene having alternative splicing
my $numAlt=0;
my @substrates=keys %records;
my $n=@substrates;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $substrate=$substrates[$i];
    my $recordsOnSubstrate=$records{$substrate};
    next unless defined $recordsOnSubstrate;
    my $nn=@$recordsOnSubstrate;
    for(my $j=0 ; $j<$nn ; ++$j)
      {
        my $record = $recordsOnSubstrate->[$j];
        my $gene = $record->[3];
        if ($gene->getNumTranscripts() > 1)
          {$gene->{alt}=1;$numAlt++}
      }
  }

################ PROCESS FASTA FILE #################
my $fastaReader;
if($opt_c)
  {$fastaReader=new CompiledFasta($fastaFile)}
else
  {$fastaReader=new FastaReader($fastaFile)}

my $fastaWriter=new FastaWriter;
my $nextId=1;
my %enumeration;
my %debugHash;
my $count;
while(1)
  {
    last unless $count<$maxN;
    my ($defline,$seq);
    if($opt_c)
      {$defline=">$opt_c"}
    else
      {($defline,$seq)=$fastaReader->nextSequenceRef()}
    last unless defined $defline;
    my $seqLen=($opt_c ? FileSize::fileSize($fastaFile) : length $$seq);
    $defline=~/>(\S+)/ || die $defline;
    my $bac="\L$1";
    my $records=$records{$bac};
    if(!defined($records))
       {
	 if($opt_c) {last}
	 next;
       }
    my $n=@$records;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	last unless $count<$maxN;
	my $record=$records->[$i];
	my ($begin,$len,$strand,$gene)=@$record;
	next if(!($gene->{alt}));

	if($begin<0) {$begin=0}
	my $end=$begin+$len;#+3;
	if($end>$seqLen) {$end=$seqLen}
	$len=$end-$begin;

	my $geneId=$gene->getId();
	if($geneId=~/(.*);/) {$geneId=$1}
	if(!defined($enumeration{$geneId})) 
	  {$enumeration{$geneId}=$nextId++}

	$debugHash{$enumeration{$geneId}}=$geneId;
	my $d=$enumeration{$geneId};

        ###### messing around $$$$$$$$$$
        my $numTranscripts = $gene->getNumTranscripts();
        #print("(round 2) numTranscripts = $numTranscripts\n");
        ###############################

	my $chunkId=$enumeration{$geneId};
	my $outfasta="$chunkId.fasta";
	my $outgff="$chunkId.gff";
        my $outcoords="$chunkId.coords.gff";
        open(OUTCOORDS,">chunks/$outcoords") || die "can't create $outcoords";
	open(OUTFASTA,">chunks/$outfasta") || die "can't create $outfasta";
	open(OUTGFF,">chunks/$outgff") || die "can't create $outgff";
        my $strand = $gene->getStrand();
        my $chr = $gene->getSubstrate();
        my $outBegin = $begin + 1; #convert to 1-based
        my $outEnd = $end; #internal coord system appears to be 0-based half-open, so end coord is the same in 1-based
        print OUTCOORDS "$chr\tarab\tchunk\t$outBegin\t$outEnd\t.\t$strand\t0\tgene_id=$geneId; chunk_id=$chunkId\n";

	my $string;
	if($opt_c)
	  {$string=$fastaReader->load($begin,$len)}
	else
	  {$string=substr($$seq,$begin,$len)}
	my $substrate=$gene->getSubstrate();
	$fastaWriter->addToFasta(">$chunkId ($substrate)",$string,
				 \*OUTFASTA);
	#$gene->setSubstrate($chunkId); do we want to do this to the transcripts? probably not

        for(my $k = 0; $k < $gene->getNumTranscripts(); $k++)
          {
            my $transcript = $gene->getIthTranscript($k);
	    # Shift coords and output GFF
	    my $numExons=$transcript->numExons();
	    my $delta=-$begin;
	    for(my $j=0 ; $j<$numExons ; ++$j)
	      {
	        my $exon=$transcript->getIthExon($j);
	        $exon->shiftCoords($delta);
	        my $gff=$exon->toGff();
	        print OUTGFF "$gff";
	      }
          }
        close(OUTGFF);
	close(OUTFASTA);
	++$count;
      }
    if($opt_c) {last}
  }

