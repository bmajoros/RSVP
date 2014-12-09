#!/usr/bin/perl

my $usage="$0 <in.gff> <out.gff> <*.fasta>";
die "$usage\n" unless @ARGV==3;
my ($infile,$outfile,$fasta)=@ARGV;

`gff-genes-to-graph.pl $infile $fasta > tmp.738271`;
`gff-graph-to-genes.pl tmp.738271 > $outfile`;




