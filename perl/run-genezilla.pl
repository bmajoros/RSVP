#!/usr/bin/perl
use strict;
use ProgramName;
use SGE;

my $EXTRA_FLAGS="";#"-O -D -S -C"; #"-D -S -C"; #"-S -C";

#my $AUTO_N="-a 1";
#my $AUTO_N="";
my $WANT_INTRONS="";#"-i";
my $WANT_UTR=0;
my $pwd=`pwd`;
chomp $pwd;
my $dashI="-i"; # this -i is for reweight-graph (one gene only)
my $dashE="";#"-e"; # reweight-graph: exon support based only on splicing
#my $DASH_Q="-q 50";
my $MAX_JOBS=1000;
my $JOB_TAG="RSVP";
my $CHUNKS_PER_QFILE=20;

my $name=ProgramName::get();
die "
$name <chunks-dir> <*.iso> <dna-weight> <rna-weight> <#predictions> <package-dir> <strict-support:0/1> <queue-size> <max-chunks>
<autoN:0/1>   (output will be in chunks-dir/*.nbest)

" unless @ARGV==10;
my ($chunks,$iso,$dna,$rna,$N,$package,$strictSupport,$queueCapacity,
    $maxChunks,$wantAutoN)=@ARGV;
if($maxChunks eq ".") { $maxChunks=100000000 }
my $DASH_Q="-q $queueCapacity";
my $sge=new SGE;
my $dashP=$strictSupport ? "-P $chunks" : "";
my $AUTO_N=($wantAutoN ? "-a 1" : "");

my @chunks;
my @files=`ls $chunks`;
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i];
  if($file=~/(\d+)\.fasta/) { push @chunks,$1 }
  last unless @chunks<$maxChunks;
}

system("rm -r q ; mkdir -p q");
system("rm -r mpirun ; mkdir -p mpirun");
system("rm -f $chunks/*.nbest");
system("rm -f $chunks/*.abinitio");
system("rm -f $chunks/*.genezilla");
system("rm -f $chunks/*.graph");
system("rm -f $chunks/*.graph2");

# #\$ -l scr_free=1G
# #\$ -q *\@ohler*,*\@igspnih-n1*,*\@igspnih-n3*,*\@igspnih-n4*,*\@igspnih-n6*,*\@igspncbc*

my $numChunks=@chunks;
for(my $j=0 ; $j<$numChunks ; ) {
  open(OUT,">q/$j.q") || die;
    print OUT "#!/bin/tcsh
#
#\$ -S /bin/tcsh -cwd
#\$ -o $pwd/mpirun/$j.mpirun -j y
#\$ -l highprio
#\$ -N $JOB_TAG$j
cd $pwd
hostname
date

";
  my $lastI;
  for(my $i=$j ; $i<$j+$CHUNKS_PER_QFILE && $i<$numChunks ; ++$i) {
    my $chunkID=$chunks[$i];
    my $dashU=$WANT_UTR ? "-u $chunks/$chunkID.pileup" : "";
    print OUT "
# CHUNK $i:

$package/c++/genezilla/genezilla $EXTRA_FLAGS -t $dashP $iso $chunks/$chunkID.fasta -o $chunks/$chunkID.graph >& $chunks/$chunkID.genezilla

#$package/c++/genezilla/genezilla $iso $chunks/$chunkID.fasta >& $chunks/$chunkID.abinitio

$package/c++/genezilla/reweight-graph $dashE $dashI $chunks/$chunkID.graph $chunks/$chunkID.pileup $chunks/$chunkID.junctions $dna $rna $chunks/$chunkID.graph2

rm $chunks/$chunkID.graph

$package/c++/genezilla/n-best $AUTO_N $WANT_INTRONS $dashU $DASH_Q $chunks/$chunkID.graph2 $N > $chunks/$chunkID.nbest
#$package/c++/genezilla/n-best $WANT_INTRONS $dashU $DASH_Q $chunks/$chunkID.graph $N > $chunks/$chunkID.nbest

#$package/perl/retain-introns.pl $chunks/$chunkID.nbest $chunks/$chunkID.pileup $chunks/$chunkID.nbest-new $package ; mv $chunks/$chunkID.nbest-new $chunks/$chunkID.nbest

rm $chunks/$chunkID.graph2

";
    $lastI=$i;
  }
  $j=$lastI+1;
  print OUT "
echo done.
date
";
}
$sge->subAllDir("q",$MAX_JOBS,$JOB_TAG);
#print "all jobs submitted.\n";




