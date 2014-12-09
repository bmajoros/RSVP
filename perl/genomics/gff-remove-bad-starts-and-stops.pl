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

# Remove genes not having canonical start/stop codons
my $n=@$transcripts;
print "read $n transcripts\n";
my %stops; $stops{"TAG"}=$stops{"TGA"}=$stops{"TAA"}=1;
my $victims=0;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
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
	undef $transcripts->[$i];
	++$victims;
      }
  }

# Write output file
open(OUT,">$outGff") || die "Can't write to $outGff\n";
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    next unless defined $transcript;
    my $gff=$transcript->toGff();
    print OUT "$gff";
  }
close(OUT);

# Report results
my $remaining=$n-$victims;
print "removed $victims transcripts due to start/stop codon anomalies\n";
print "remaining transcripts: $remaining\n";





