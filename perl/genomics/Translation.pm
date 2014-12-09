package Translation;
use strict;

initCodonMap();
initComplementMap();

######################################################################
#
# Translation.pm
#
# bmajoros@tigr.org 7/8/2001
#
# Routines for translation and reverse complementing of nucleotide
# sequences.
#
# Attributes:
#
# Methods:
#    $aaSeq=Translation::translate(\$nucSeq);
#    $revSeq=Translation::reverseComplement(\$nucSeq);
#
# Private methods:
#    initCodonMap();
#    initComplementMap();
#    $nucleotide=$complement($nucleotide);
######################################################################


#---------------------------------------------------------------------
#                           PUBLIC METHODS
#---------------------------------------------------------------------
#    $aaSeq=Translation::translate(\$nucSeq);
sub translate
{
  my ($transcript)=@_;

  my $translation;

  my $i;
  my $len=length $$transcript;
  for($i=0 ; $i+2<$len ; $i+=3)
    {
      my $three=substr($$transcript,$i,3);
      my $residue=$::codon{$three};
      $residue="X" unless defined $residue;
      $translation.=$residue;
      #last if $residue eq "*";
    }

  return $translation;
}
#---------------------------------------------------------------------
#    $revSeq=Translation::reverseComplement($nucSeq);
sub reverseComplement
{
  my ($seq)=@_;
  my $i;
  my $len=length $$seq;
  my $buffer="";
  for($i=0 ; $i<$len ; ++$i)
    {
      $buffer.=complement(substr($$seq,$len-1-$i,1));
    }
  return $buffer;
}
#---------------------------------------------------------------------










#---------------------------------------------------------------------
#                         PRIVATE METHODS
#---------------------------------------------------------------------
sub initComplementMap
{
  $::complementMap{'A'}= 'T';
  $::complementMap{'T'}= 'A';
  $::complementMap{'G'}= 'C';
  $::complementMap{'C'}= 'G';
  $::complementMap{'R'}= 'Y';
  $::complementMap{'Y'}= 'R';

  $::complementMap{'N'}= 'N'; ### <------DEBUGGING!
}
#---------------------------------------------------------------------
sub complement
{
  my ($c)=@_;
  if(!defined($::complementMap{$c})) 
    {
      $::complementMap{$c}="N";
      #die "base $c has no complement"
    }
  return $::complementMap{$c};
}
#---------------------------------------------------------------------
sub initCodonMap
{
  $::codon{"GTT"}='V';
  $::codon{"GTC"}='V';
  $::codon{"GTA"}='V';
  $::codon{"GTG"}='V';
  $::codon{"GTR"}='V';
  $::codon{"GTY"}='V';
  $::codon{"GTN"}='V';
  $::codon{"GCT"}='A';
  $::codon{"GCC"}='A';
  $::codon{"GCA"}='A';
  $::codon{"GCG"}='A';
  $::codon{"GCR"}='A';
  $::codon{"GCY"}='A';
  $::codon{"GCN"}='A';
  $::codon{"GAT"}='D';
  $::codon{"GAC"}='D';
  $::codon{"GAY"}='D';
  $::codon{"GAA"}='E';
  $::codon{"GAG"}='E';
  $::codon{"GAR"}='E';
  $::codon{"GGT"}='G';
  $::codon{"GGC"}='G';
  $::codon{"GGA"}='G';
  $::codon{"GGG"}='G';
  $::codon{"GGR"}='G';
  $::codon{"GGY"}='G';
  $::codon{"GGN"}='G';
  $::codon{"TTT"}='F';
  $::codon{"TTC"}='F';
  $::codon{"TTY"}='F';
  $::codon{"TTA"}='L';
  $::codon{"TTG"}='L';
  $::codon{"TTR"}='L';
  $::codon{"CTT"}='L';
  $::codon{"CTC"}='L';
  $::codon{"CTA"}='L';
  $::codon{"CTG"}='L';
  $::codon{"CTN"}='L';
  $::codon{"YTR"}='L';
  $::codon{"TCT"}='S';
  $::codon{"TCC"}='S';
  $::codon{"TCA"}='S';
  $::codon{"TCG"}='S';
  $::codon{"TCN"}='S';
  $::codon{"AGT"}='S';
  $::codon{"AGC"}='S';
  $::codon{"AGY"}='S';
  $::codon{"TAT"}='Y';
  $::codon{"TAC"}='Y';
  $::codon{"TAY"}='Y';
  $::codon{"TAA"}='*';
  $::codon{"TAG"}='*';
  $::codon{"TAR"}='*';
  $::codon{"TAG"}='*';
  $::codon{"TGT"}='C';
  $::codon{"TGC"}='C';
  $::codon{"TGY"}='C';
  $::codon{"TGA"}='*';
  $::codon{"TGG"}='W';
  $::codon{"CCT"}='P';
  $::codon{"CCC"}='P';
  $::codon{"CCA"}='P';
  $::codon{"CCG"}='P';
  $::codon{"CCR"}='P';
  $::codon{"CCY"}='P';
  $::codon{"CCN"}='P';
  $::codon{"CAT"}='H';
  $::codon{"CAC"}='H';
  $::codon{"CAY"}='H';
  $::codon{"CAA"}='Q';
  $::codon{"CAG"}='Q';
  $::codon{"CAR"}='Q';
  $::codon{"CGT"}='R';
  $::codon{"CGC"}='R';
  $::codon{"CGA"}='R';
  $::codon{"CGG"}='R';
  $::codon{"CGN"}='R';
  $::codon{"ATT"}='I';
  $::codon{"ATC"}='I';
  $::codon{"ATA"}='I';
  $::codon{"ATH"}='I';
  $::codon{"ATG"}='M';
  $::codon{"ACT"}='T';
  $::codon{"ACC"}='T';
  $::codon{"ACA"}='T';
  $::codon{"ACG"}='T';
  $::codon{"ACN"}='T';
  $::codon{"AAT"}='N';
  $::codon{"AAC"}='N';
  $::codon{"AAY"}='N';
  $::codon{"AAA"}='K';
  $::codon{"AAG"}='K';
  $::codon{"AAR"}='K';
  $::codon{"AGG"}='R';
  $::codon{"AGA"}='R';
  $::codon{"AGY"}='R';
}
#---------------------------------------------------------------------

1;

