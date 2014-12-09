#!/usr/bin/perl
use strict;
use GenbankParser;
use FastaWriter;

$0=~/([^\/]+$)/;
my $program=$1;
my $usage="$program <infile> <gff-outfile> <fasta-outfile>\n";
die $usage unless @ARGV==3;
my ($infile,$gffFile,$fastaFile)=@ARGV;

my $fastaWriter=new FastaWriter;
open(GFF,">$gffFile") || die "Can't create file: $gffFile\n";
open(FASTA,">$fastaFile") || die "Can't create file: $fastaFile\n";
my $transgrp=1;
my $parser=new GenbankParser($infile);
while(1)
  {
    my $entry=$parser->nextEntry();
    last unless $entry;

    my $defline=">$transgrp";
    my $organism=$entry->findUnique("ORGANISM");
      # || die "No ORGANISM field in Genbank entry!\n";
    if(!defined($organism)) {$organism="unknown_organism"}
    my @words=split/\s+/,$organism;
    $organism="$words[0]_$words[1]";
    my $seq=$entry->getSubstrate();

    my $geneId=$transgrp;
    my $locus=$entry->findUnique("LOCUS");
    if($locus) {$geneId=$locus}

    my $features=$entry->findUnique("FEATURES");
    if($features)
      {
	my $geneElem=$features->findUnique("gene");
	if($geneElem)
	  {if($geneElem=~/\/gene=\"([^\"]+)\"/) {$geneId=$1}}

	my $fivePrimeUTR=$features->findUnique("5'UTR");
	if($fivePrimeUTR=~/^\s*(complement\()?(<)?(\d+)\.\.(>)?(\d+)/)
	  {
	    my ($complement,$less,$begin,$greater,$end)=($1,$2,$3,$4,$5);
	    my $strand=($begin<$end ? "+" : "-");
	    if($complement) {$strand=($strand eq "+" ? "-" : "+")}
	    if($strand eq "-") {($begin,$end)=($end,$begin)}
	    print GFF "$transgrp\t$organism\tfive-prime-UTR\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	  }
	my $threePrimeUTR=$features->findUnique("3'UTR");
	if($threePrimeUTR=~/^\s*(complement\()?(<)?(\d+)\.\.(>)?(\d+)/)
	  {
	    my ($complement,$less,$begin,$greater,$end)=($1,$2,$3,$4,$5);
	    my $strand=($begin<$end ? "+" : "-");
	    if($complement) {$strand=($strand eq "+" ? "-" : "+")}
	    if($strand eq "-") {($begin,$end)=($end,$begin)}
	    print GFF "$transgrp\t$organism\tthree-prime-UTR\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	  }
	my $signalPeptide=$features->findUnique("sig_peptide");
	if($signalPeptide=~/^\s*(complement\()?(<)?(\d+)\.\.(>)?(\d+)/)
	  {
	    my ($complement,$less,$begin,$greater,$end)=($1,$2,$3,$4,$5);
	    my $strand=($begin<$end ? "+" : "-");
	    if($complement) {$strand=($strand eq "+" ? "-" : "+")}
	    if($strand eq "-") {($begin,$end)=($end,$begin)}
	    print GFF "$transgrp\t$organism\tsignal_peptide\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	  }
	elsif($signalPeptide=~/^\s*join\(([\d\.,]+)\)/)
	  {
	    my @fields=split/,/,$1;
	    foreach my $field (@fields)
	      {
		$field=~/(\d+)\.\.(\d+)/ || die $field;
		my ($begin,$end)=($1,$2);
		my $strand=($begin<$end ? "+" : "-");
		if($strand eq "-") {($begin,$end)=($end,$begin)}
		print GFF "$transgrp\t$organism\tsignal_peptide\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	      }
	  }
	my $tata=$features->findUnique("TATA_signal");
	if($tata=~/^\s*(complement\()?(<)?(\d+)\.\.(>)?(\d+)/)
	  {
	    my ($complement,$less,$begin,$greater,$end)=($1,$2,$3,$4,$5);
	    my $strand=($begin<$end ? "+" : "-");
	    if($complement) {$strand=($strand eq "+" ? "-" : "+")}
	    if($strand eq "-") {($begin,$end)=($end,$begin)}
	    print GFF "$transgrp\t$organism\tTATA\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	  }
	my $polyAs=$features->findPairs("polyA_signal");
	foreach my $pair (@$polyAs)
	  {
	    my $polyA=$pair->[1];
	    if($polyA=~/^\s*(complement\()?(<)?(\d+)\.\.(>)?(\d+)/)
	      {
		my ($complement,$less,$begin,$greater,$end)=
		  ($1,$2,$3,$4,$5);
		my $strand=($begin<$end ? "+" : "-");
		if($complement) {$strand=($strand eq "+" ? "-" : "+")}
		if($strand eq "-") {($begin,$end)=($end,$begin)}
		print GFF "$transgrp\t$organism\tpolyA\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	      }
	  }
      }

    my $exons=$entry->getCDS();
    unless($exons) {++$transgrp;next}
    #die "No CDS in Genbank entry!\n" unless $exons;
    my $numExons=@$exons;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my ($begin,$end)=@$exon;
	my $strand=($begin<$end ? "+" : "-");
	if($strand eq "-") {($begin,$end)=($end,$begin)}
	print GFF "$transgrp\t$organism\tCDS\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	$fastaWriter->addToFasta($defline,$seq,\*FASTA);
      }

    $exons=$entry->getExons();
    unless($exons) {++$transgrp;next}
    $numExons=@$exons;
    for(my $i=0 ; $i<$numExons ; ++$i)
      {
	my $exon=$exons->[$i];
	my ($begin,$end)=@$exon;
	my $strand=($begin<$end ? "+" : "-");
	if($strand eq "-") {($begin,$end)=($end,$begin)}
	print GFF "$transgrp\t$organism\texon\t$begin\t$end\t.\t$strand\t.\ttranscript_id=\"$transgrp\"; gene_id=\"$geneId\"\n";
	$fastaWriter->addToFasta($defline,$seq,\*FASTA);
      }
    ++$transgrp;
  }
close(FASTA);
close(GFF);




