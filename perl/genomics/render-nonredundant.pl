#!/usr/bin/perl
use strict;
use BlastParser;
use ProgramName;
$|=1;

my $name=ProgramName::get();
my $usage="$name <blast-output> <max-expect> <max-\%aligned> <max-\%identity>

example:
  $0 blast.out 0.001 0.95 0.90

";
die "$usage\n" unless @ARGV==4;
my ($blastFile,$maxExpect,$maxAligned,$maxIdentity)=@ARGV;

my %allQueryIds;
my $parser=new BlastParser;
my $hits=$parser->parse($blastFile,$maxExpect);
my $nHits=@$hits;
my %signatures;
print STDERR "$nHits hits with E<=$maxExpect\n";
for(my $i=0 ; $i<$nHits ; ++$i)
  {
    my $hitPair=$hits->[$i];
    my ($queryId,$match)=@$hitPair;
    $allQueryIds{$queryId}=1;
    my $matchId=$match->{id};
    next if $queryId eq $matchId;
    my $signature="$queryId\=$matchId";
    push @{$signatures{$signature}},$hitPair;
  }
my @allQueryIds=keys %allQueryIds;
my $numQueries=@allQueryIds;
print STDERR "$numQueries query sequences\n";

my @signatures=keys %signatures;
my $n=@signatures;
for(my $i=0 ; $i<$n ; ++$i)
  {
    my $signature=$signatures[$i];
    $signature=~/(.+)=(.+)/ || die "can't parse: $signature";
    my ($queryId,$matchId)=($1,$2);
    my $HSPs=$signatures{$signature};
    my $nHSPs=@$HSPs;

    my ($totalIdentities,$totalAlignedLength,$queryLength);
    for(my $j=0 ; $j<$nHSPs ; ++$j)
      {
	my $hitPair=$HSPs->[$j];
	my $matchInfo=$hitPair->[1];
	$queryLength=$matchInfo->{queryLength};
	$totalAlignedLength+=$matchInfo->{alignLength};
	$totalIdentities+=$matchInfo->{identities};
      }
    my $percentAligned=$totalAlignedLength/$queryLength;
    my $percentIdentity=$totalIdentities/$totalAlignedLength;

    my $pctAligned=int($percentAligned*100+5/9);
    my $pctIdentity=int($percentIdentity*100+5/9);
    #print STDERR "$queryId=$matchId, $nHSPs HSPs, $pctAligned\% aligned, $pctIdentity\% identity\n";
    if($percentAligned>=$maxAligned && $percentIdentity>=$maxIdentity)
      {
	print STDERR "SIGNIFICANT MATCH: $queryId=$matchId ($percentIdentity over $percentAligned)\n";
	if($allQueryIds{$queryId} && $allQueryIds{$matchId}) 
	  {undef $allQueryIds{$queryId};if(defined($allQueryIds{$queryId})){die}}
      }
  }
my @survivors;#=keys %allQueryIds;
for(my $i=0 ; $i<$numQueries ; ++$i)
  {
    my $id=$allQueryIds[$i];
    if(defined($allQueryIds{$id})) {push @survivors,$id}
  }
my $numSurvivors=@survivors;
my $percentSurviving=int(100*$numSurvivors/$numQueries+5/9);
print STDERR "$percentSurviving\% remain ($numSurvivors/$numQueries)\n";
for(my $i=0 ; $i<$numSurvivors ; ++$i)
  {print "$survivors[$i]\n"}







