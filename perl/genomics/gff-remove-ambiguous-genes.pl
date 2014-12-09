#!/usr/bin/perl
use strict;
use GffTranscriptReader;
use CompiledFasta;
use Translation;

# Process command line and load GFF file
$0=~/([^\/]+)\s*$/;
my $usage="$1 <in.gff> <compiled-fasta-file> <out.gff>";
die "$usage\n" unless @ARGV==3;
my ($inGff,$fasta,$outGff)=@ARGV;
my $fastaReader=new CompiledFasta($fasta);
my $gffReader=new GffTranscriptReader;
my $transcripts=$gffReader->loadGFF($inGff);

# Remove genes having alternative transcripts
my %transcriptsPerGene;
my $numTranscripts=@$transcripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $geneId=$transcript->getGeneId();
    ++$transcriptsPerGene{$geneId};
  }
my @goodTranscripts;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $geneId=$transcript->getGeneId();
    if($transcriptsPerGene{$geneId}==1) {push @goodTranscripts,$transcript}
  }
undef $transcripts;
my $n=@goodTranscripts;
print "#transcripts reduced from $numTranscripts to $n due to alternative splicing\n";

# Remove genes not having canonical start/stop codons
my $n=@goodTranscripts;
print "read $n transcripts\n";
my %stops; $stops{"TAG"}=$stops{"TGA"}=$stops{"TAA"}=1;
my $victims=0;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$goodTranscripts[$i];
    my $geneId=$transcript->getGeneId();
    my $strand=$transcript->getStrand();
    my $begin=$transcript->getBegin();
    my $end=$transcript->getEnd();
    my $transSeq=$transcript->loadTranscriptSeqCF($fastaReader);
    my $startCodon=substr($transSeq,0,3);
    my $stopCodon=substr($transSeq,length($transSeq)-3,3);
    $startCodon="\U$startCodon";
    $stopCodon="\U$stopCodon";
    if($startCodon ne "ATG" || !defined($stops{$stopCodon})) 
      {
	undef $goodTranscripts[$i];
	++$victims;
      }
  }

# Write output file
open(OUT,">$outGff") || die "Can't write to $outGff\n";
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$goodTranscripts[$i];
    next unless defined $transcript;
    my $gff=$transcript->toGff();
    print OUT "$gff";
  }
close(OUT);

# Report results
my $remaining=$n-$victims;
print "removed a further $victims transcripts due to start/stop codon anomalies\n";
print "remaining transcripts: $remaining\n";





