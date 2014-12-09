#!/usr/bin/perl
use strict;
use FastaReader;
use FastaWriter;
use GffTranscriptReader;

# Process command line
my $usage="$0 <*.gff> <*.fasta> <max-size> <output-filestem>";
die "$usage\n" unless @ARGV==4;
my ($inGff,$inFasta,$maxSize,$filestem)=@ARGV;

# Read transcripts
my $gffReader=new GffTranscriptReader();
my $transcripts=$gffReader->loadGFF($inGff);
my $numTranscripts=@$transcripts;

# Enter the transcripts into a hash table
my %transcriptsByContig;
for(my $i=0 ; $i<$numTranscripts ; ++$i)
  {
    my $transcript=$transcripts->[$i];
    my $id=$transcript->getSubstrate();
    push @{$transcriptsByContig{$id}},$transcript;
  }

# Create output files and writers
my $writer=new FastaWriter();
my $outGff="$filestem.gff";
my $outFasta="$filestem.fasta";
open(OUTGFF,">$outGff") || die "Can't create file $outGff\n";
open(OUTFASTA,">$outFasta") || die "Can't create file $outFasta\n";

# Process all input sequences
my ($totalTranscripts,$lostTranscripts,$totalChunks);
my $fastaReader=new FastaReader($inFasta);
while(1)
  {
    # Get the next sequence and all the transcripts on that substrate
    my ($def,$seq)=$fastaReader->nextSequence();
    last unless defined $def;
    $def=~/^\s*>\s*(\S+)(.*)/ || die "can't parse defline: $def\n";
    my ($id,$rest)=($1,$2);
    my $transcripts=$transcriptsByContig{$id};
    my $numTranscripts=@$transcripts;
    $totalTranscripts+=$numTranscripts;

    # Repeatedly chop pieces off the front, until no part of the sequence
    # is too long
    my $len=length $seq;
    my $nextId=1;
    while($len>$maxSize)
      {
	# Chop off the front and output it as a separate "contig"
	my $piece=substr($seq,0,$maxSize);
	$len-=$maxSize;
	$seq=substr($seq,$maxSize,$len);
	my $newId="${nextId}_$id";
	++$nextId;
	my $defline=">$newId $rest";
	$writer->addToFasta($defline,$piece,\*OUTFASTA);
	++$totalChunks;
	for(my $i=0 ; $i<$numTranscripts ; ++$i)
	  {
	    my $transcript=$transcripts->[$i];
	    next unless defined $transcript;
	    if($transcript->isContainedWithin(0,$maxSize))
	      {
		$transcript->setSubstrate($newId);
		my $gff=$transcript->toGff();
		print OUTGFF "$gff";
		undef $transcripts->[$i];
	      }
	    else
	      {
		# Adjust the coordinates of the transcripts on the remaining
		# portion of the contig
		$transcript->shiftCoords(-$maxSize);
	      }
	  }
      }

    # Output the remaining part of the sequence & its transcripts
    if($len>0)
      {
	my $newId="${nextId}_$id";
	my $defline=">$newId $rest";
	$writer->addToFasta($defline,$seq,\*OUTFASTA);
	++$totalChunks;
	for(my $i=0 ; $i<$numTranscripts ; ++$i)
	  {
	    my $transcript=$transcripts->[$i];
	    next unless defined $transcript;
	    if($transcript->isContainedWithin(0,$maxSize))
	      {
		$transcript->setSubstrate($newId);
		my $gff=$transcript->toGff();
		print OUTGFF "$gff";
		#undef $transcripts->[$i];
	      }
	    else {++$lostTranscripts}
	  }
      }
  }

close(OUTFASTA);
close(OUTGFF);

print "$totalChunks total sequences written\n";
print "$lostTranscripts transcripts (of $totalTranscripts) lost due to straddling a split\n";


