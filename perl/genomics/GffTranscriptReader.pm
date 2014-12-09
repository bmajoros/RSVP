package GffTranscriptReader;
use strict;
use Exon;
use Transcript;
use Gene;

######################################################################
#
# GffTranscriptReader.pm
#
# bmajoros@tigr.org 7/8/2002
#
# Returns a list of Transcripts.  For each transcript, the Exons
# will be sorted according to order of translation, so that
# the exon containing the start codon will come before the exon
# containing the stop codon.  This means that for the minus strand,
# exon begin coordinates will be decreasing.  However, an individual
# exon's begin coordinate is always less than its end coordinate.
# The transcripts themselves are sorted along the chromosome left-to-
# right.  Note that although the GFF coordinates are 1-based/base-based
# (1/B), internally all coordinates are stored as 0-based/space-based
# (0/S).  This conversion is handled automatically.
#
# Attributes:
#   shouldSortTranscripts
# Methods:
#   $reader=new GffTranscriptReader();
#   $reader->setStopCodons({TAG=>1,TAA=>1,TGA=>1});
#   $transcriptArray=$reader->loadGFF($filename);
#   $geneList=$reader->loadGenes($filename);
#   $hashTable=$reader->hashBySubstrate($filename);
#   $hashTable=$reader->loadTranscriptIdHash($filename);
#   $hashTable=$reader->loadGeneIdHash($filename);
#   $reader->doNotSortTranscripts();
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
sub new
{
  my ($class)=@_;

  my $self=
    {
     shouldSortTranscripts=>1,
     stopCodons=>{TAG=>1,TAA=>1,TGA=>1},
    };
  bless $self,$class;

  return $self;
}
#---------------------------------------------------------------------
#   $geneList=$reader->loadGenes($filename);
sub loadGenes
  {
    my ($self,$filename)=@_;
    my $transcripts=$self->loadGFF($filename);
    my $n=@$transcripts;
    my %genes;
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $transcript=$transcripts->[$i];
	my $gene=$transcript->getGene();
	$genes{$gene}=$gene;
      }
    my @genes=values %genes;
    @genes=sort {$a->getBegin() <=> $b->getBegin()} @genes;
    return \@genes;
  }
#---------------------------------------------------------------------
# $transcriptArray=$reader->loadGFF($filename);
#
sub loadGFF
  {
    my ($self,$gffFilename)=@_;
    my (%transcripts,%genes,$geneId);
    my $readOrder=1;
    if($gffFilename=~/\.gz$/)
      {open(GFF,"cat $gffFilename|gunzip|") || die $gffFilename}
    else
      {open(GFF,$gffFilename) || die $gffFilename}
    while(<GFF>)
      {
	next unless $_=~/\S+/;
	next if $_=~/^\s*\#/;
	my @fields=split/\s+/,$_;
	if($fields[2]=~/exon/ || $fields[2]=~/CDS/)
	  {
	    my $exonBegin=$fields[3]-1;
	    my $exonEnd=$fields[4];
	    my $exonScore=$fields[5];
	    my $strand=$fields[6];
	    my $frame=$fields[7];
	    my $transcriptId;#=$fields[8];
	    if($_=~/transgrp[:=]\s*(\S+)/) {$transcriptId=$1}
	    elsif($_=~/transcript_id[:=]?\s*\"?([^\s\"]+)\"?/) 
	      {$transcriptId=$1}
	    elsif($_=~/Parent=([^;,\s]+)/) {$transcriptId=$1}
	    my $geneId;#=$fields[8];
	    if(/genegrp=(\S+)/) {$geneId=$1}
	    elsif(/gene_id[:=]?\s*\"?([^\s\"]+)\"?/) {$geneId=$1}
	    if(!defined($transcriptId)) {$transcriptId=$geneId}
	    if(!defined($geneId)) {$geneId=$transcriptId}
	    chop $transcriptId if $transcriptId=~/;$/;
	    chop $geneId if $geneId=~/;$/;
	    my $extra;
	    for(my $i=8;$i<@fields;++$i){$extra.=$fields[$i]." "}
	    if($exonBegin>$exonEnd)
	      {($exonBegin,$exonEnd)=($exonEnd,$exonBegin)}
	    my $transcript=$transcripts{$transcriptId};
	    if(!defined $transcript)
	      {
		$transcripts{$transcriptId}=$transcript=
		  new Transcript($transcriptId,$strand);
		$transcript->setStopCodons($self->{stopCodons});
		$transcript->{readOrder}=$readOrder++;
		$transcript->{substrate}=$fields[0];
		#if($transcript->{substrate}=~/([^\.]+)\./)
		#  {$transcript->{substrate}=$1}
		$transcript->{source}=$fields[1];
	      }
	    $transcript->{geneId}=$geneId;
	    $transcript->{extraFields}=$extra;
	    my $gene=$genes{$geneId};
	    if(!defined $gene) 
	      {$genes{$geneId}=$gene=new Gene(); $gene->setId($geneId)}
	    $transcript->setGene($gene);
	    my $exon=new Exon($exonBegin,$exonEnd,$transcript);
	    if(!$transcript->exonOverlapsExon($exon))
	      {
		$exon->{frame}=$frame;
		$exon->{score}=$exonScore;
		$exon->{type}=$fields[2];
		
		push @{$transcript->{exons}},$exon; # OK -- we sort later
		#$transcript->addExon($exon);       # <---too slow
	      }
	    if($transcript->numExons()==1) # first exon seen...
	      {$gene->addTranscript($transcript)}
	  }
	elsif($_=~/translation\s+(\d+)\s+(\d+).*([\-\+]).*transgrp=(\S+)/ ||
	      $_=~/translation\s+(\d+)\s+(\d+).*([\-\+]).*transcript_id:\s*\"?([^\s\"]+)\"?/)
	  {
	    my ($startCodon,$stopCodon,$strand,$transcriptId)=
	      ($1-1,$2,$3,$4,);
	    $_=~/genegrp=(\S+)|gene_id:?\s*\"[\s\"]+\"/ || die $_;
	    my $geneId=$1;
	    my $transcript=$transcripts{$transcriptId};
	    if($strand eq "-") {$startCodon=$stopCodon}
	    $transcript->{startCodon}=$startCodon;
	    $transcript->{startCodonAbsolute}=$startCodon;
	    $transcript->{geneId}=$geneId;
	  }
	elsif($fields[2]=~/start-codon/)
	  {
	    my $startCodonBegin=$fields[3]-1;
	    my $startCodonEnd=$fields[4];
	    my $strand=$fields[6];
	    my $transcriptId=$fields[8];
	    if($_=~/transgrp=(\S+)/) {$transcriptId=$1}
	    elsif($_=~/transcript_id:\s*\"?([^\s\"]+)\"?/) 
	      {$transcriptId=$1}
	    my $geneId=$fields[8];
	    if(/genegrp=(\S+)/) {$geneId=$1}
	    elsif(/gene_id:?\s*\"?([^\s\"]+)\"?/) {$geneId=$1}
	    if($startCodonEnd ne "." &&
	       $startCodonBegin>$startCodonEnd)
	      {($startCodonBegin,$startCodonEnd)=
		 ($startCodonEnd,$startCodonBegin)}
	    my $transcript=$transcripts{$transcriptId};
	    if(!defined $transcript)
	      {
		$transcripts{$transcriptId}=$transcript=
		  new Transcript($transcriptId,$strand);
		$transcript->setStopCodons($self->{stopCodons});
		$transcript->{readOrder}=$readOrder++;
		$transcript->{substrate}=$fields[0];
		$transcript->{source}=$fields[1];
	      }
	    $transcript->{geneId}=$geneId;
	    $transcript->{startCodon}=$startCodonBegin;
	    $transcript->{startCodonAbsolute}=$startCodonBegin;
	  }
      }
    close(GFF);
    my $transcripts=[];
    @$transcripts=values %transcripts;
    adjustStartCodons($transcripts);
    computeFrames($transcripts);
    undef %transcripts;
    if($self->{shouldSortTranscripts})
      {
	@$transcripts=sort
	  {
	    my $cmp=0+($a->{substrate} cmp $b->{substrate});
	    if($cmp==0) {$cmp=($a->{begin}<=>$b->{begin})}
	    $cmp;
	  } @$transcripts;
      }
    else 
      {@$transcripts=sort {$a->{readOrder}<=>$b->{readOrder}} @$transcripts}
    return $transcripts;
  }
#--------------------------------------------------------------------------
#   $reader->doNotSortTranscripts();
sub doNotSortTranscripts
  {
    my ($self)=@_;
    $self->{shouldSortTranscripts}=0;
  }
#--------------------------------------------------------------------------
#   $hashTable=$reader->loadTranscriptIdHash($filename);
sub loadTranscriptIdHash
  {
    my ($self,$filename)=@_;
    my $transcriptArray=$self->loadGFF($filename);
    my $n=@$transcriptArray;
    my $hash={};
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $transcript=$transcriptArray->[$i];
	my $id=$transcript->getID();
	$hash->{$id}=$transcript;
      }
    return $hash;
  }
#--------------------------------------------------------------------------
#   $hashTable=$reader->loadGeneIdHash($filename);
sub loadGeneIdHash
  {
    my ($self,$filename)=@_;
    my $transcriptArray=$self->loadGFF($filename);
    my $n=@$transcriptArray;
    my $hash={};
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $transcript=$transcriptArray->[$i];
	my $id=$transcript->getGeneId();
	push @{$hash->{$id}},$transcript;
      }
    return $hash;
  }
#--------------------------------------------------------------------------
#   $hashTable=$reader->hashBySubstrate($filename);
sub hashBySubstrate
  {
    my ($self,$filename)=@_;
    my $transcriptArray=$self->loadGFF($filename);
    my $n=@$transcriptArray;
    my $hash={};
    for(my $i=0 ; $i<$n ; ++$i)
      {
	my $transcript=$transcriptArray->[$i];
	my $id=$transcript->getSubstrate();
	push @{$hash->{$id}},$transcript;
      }
    return $hash;
  }
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------






#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
sub computeFrames
{
    my ($transcripts)=@_;
    foreach my $transcript (@$transcripts)
    {
      if(!$transcript->areExonTypesSet()) {$transcript->setExonTypes()}
      my $strand=$transcript->getStrand();
      my $exons=$transcript->{exons};
      my $n=@$exons;
      my $frame=0; # fine for both strands
      for(my $i=0 ; $i<$n ; ++$i)
	{
	  my $exon=$exons->[$i];
	  $exon->{frame}=$frame;
	  my $length=$exon->getLength();
	  $frame=($frame+$length)%3; # this is fine, on both strands
	}
    }
}
#---------------------------------------------------------------------
sub adjustStartCodons
{
    my ($transcripts)=@_;
    foreach my $transcript (@$transcripts)
    {
      $transcript->sortExons();
      $transcript->adjustOrders();
	my $strand=$transcript->{strand};
	my $startCodon;
	my $totalIntronSize=0;
	if($strand eq "+")
	{
	    my @exons=
		sort {$a->{begin} <=> $b->{begin}} @{$transcript->{exons}};
	    $transcript->{exons}=\@exons;
	    my $numExons=@exons;
	    next unless $numExons>0;
	    $transcript->{begin}=$exons[0]->{begin};
	    $transcript->{end}=$exons[$numExons-1]->{end};
	    if(defined($transcript->{startCodon}))
	      {
		$startCodon=$transcript->{startCodon}-$transcript->{begin};
	      }
	    else
	      {
		$startCodon=0;
		$transcript->{startCodon}=$transcript->{begin};
		$transcript->{startCodonAbsolute}=$transcript->{begin};
	      }
	    for(my $i=0 ; $i<$numExons ; ++$i)
	      {
		my $exon=$exons[$i];
		$exon->{order}=$i;
	      }
	    for(my $i=0 ; $i<$numExons ; ++$i)
	    {
		my $exon=$exons[$i];
		if($i>0)
		{
		    my $prevExon=$exons[$i-1];
		    my $intronSize=$exon->{begin}-$prevExon->{end};
		    $totalIntronSize+=$intronSize;
		}
		if(defined($transcript->{startCodon}) && 
		   exonContainsPoint($exon,$transcript->{startCodon})) 
		  {last}
	    }
	}
	else # $strand eq "-"
	{
            my @exons=
                sort {$b->{begin} <=> $a->{begin}} @{$transcript->{exons}};
            $transcript->{exons}=\@exons;
            my $numExons=@exons;
	    next unless $numExons>0;
            $transcript->{end}=$exons[0]->{end};
            $transcript->{begin}=$exons[$numExons-1]->{begin};

	    if(defined $transcript->{startCodon})
	      {$startCodon=$transcript->{end}-$transcript->{startCodon}}
	    else
	      {
		$startCodon=0;
		$transcript->{startCodon}=$transcript->{end};###
		$transcript->{startCodonAbsolute}=$transcript->{end};###
	      }
            for(my $i=0 ; $i<$numExons ; ++$i)
            {
                my $exon=$exons[$i];
		$exon->{order}=$i;
                if($i>0)
                {
                    my $prevExon=$exons[$i-1];
                    my $intronSize=$prevExon->{begin}-$exon->{end};
                    $totalIntronSize+=$intronSize;
                }
                if(defined($transcript->{startCodon}) &&
		   exonContainsPoint($exon,$transcript->{startCodon})) 
		  {last}
            }
	}
	if(defined($startCodon))
	  {
	    $startCodon-=$totalIntronSize;
	    $transcript->{startCodon}=$startCodon;
	  }
    }
}
#--------------------------------------------------------------------------
sub exonContainsPoint
{
    my ($exon,$point)=@_;
    return $point>=$exon->{begin} &&
      $point<=$exon->{end}; # this '<=' is necessary for the minus strand!
                            # do not change it back to '<' !!!
}
#--------------------------------------------------------------------------
#   $reader->setStopCodons({TAG=>1,TAA=>1,TGA=>1});
sub setStopCodons
  {
    my ($self,$stopCodons)=@_;
    $self->{stopCodons}=$stopCodons;
  }
#--------------------------------------------------------------------------

1;

