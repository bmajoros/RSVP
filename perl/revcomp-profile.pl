#!/usr/bin/perl
#Reverse complement the read and splice profiles produced by get-read-profiles,
#so that they can be used in forward strand gene prediction.

use warnings;
use strict;

die "revcomp-profile.pl <profile-file>" unless @ARGV==1;

my $file = shift;
open(OLD_PROFILE, $file) or die "profile not found";
my @oldProfile = <OLD_PROFILE>;
close(OLD_PROFILE);

my $header = $oldProfile[0];
chomp($header);
my @headerSplit = split(/\s+/, $header);
my ($start, $end) = ($headerSplit[2], $headerSplit[3]);
my $intEnd = $end - $start; 

my @newProfile;
my $newHeader = "$header REVCOMP\n";
$newProfile[0] = $newHeader;
for (my $i = 1; $i < @oldProfile; $i++)
{
    my $count = $oldProfile[$i]; #includes newline
    #convert line number to original internal coordinate
    my $oldCoord = $i - 1;
    #revcomp internal coordinate (internal end - original)
    my $newCoord = $intEnd - $oldCoord;
    #use revcomped coordinate (+1) as index in array
    my $newIndex = $newCoord + 1;
    $newProfile[$newIndex] = $count;
}

if (@newProfile != @oldProfile) {print "error--new profile has different length than old profile\n";}

print(@newProfile);
