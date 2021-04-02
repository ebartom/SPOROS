#!/usr/bin/perl -w
use strict;
use warnings;

# This script assumes a UMI with 4N's before the sequenced species and 2N's after.

my $countfile = $ARGV[0];
my $beforeN = 4;
my $afterN = 2;
my %counts;
my %reads;
my %umis;
open(IN,$countfile);
while (<IN>){
    my $stuff = 0;
    my $count = 0;
    my $read = "";
    ($stuff,$count,$read) = split(/\s+/,$_);
    if ($stuff ne ""){ $read = $count; $count = $stuff;}
    if ((length($read) > ($beforeN + $afterN)) && ($read !~ /N/)){
	my $before = substr($read,0,$beforeN);
	my $after = substr($read,-$afterN,$afterN);
	my $RNAlength = length($read) - $beforeN - $afterN;
	my $RNA = substr($read,$beforeN,$RNAlength);
	my $UMI = $before.$after;
	#	print STDERR "$read\n$before $RNA $after\t$count\n";
	if (exists ($reads{$RNA})){
	    $reads{$RNA} += 1;
	} else {$reads{$RNA} = 1;}
	if ($UMI !~ m/N/){
	    $umis{$UMI} = "found";
	    $counts{$RNA}{$UMI} = $count;
	}
    }
}
my ($umidfile,$stdfile) = "";
if ($countfile  =~ /^(.*)\.justReads.uniqCounts.txt/){
    $umidfile = "$1.justReads.uniqCounts.UMId.txt";
    $stdfile = "$1.justReads.uniqCounts.notUMId.txt";
}
 
open(UMID,">$umidfile");
open(STD,">$stdfile");   
foreach my $RNA (keys(%reads)){
    print UMID "$RNA\t$reads{$RNA}\n";
    my $sumcount = 0;
    foreach my $UMI (keys(%umis)){
	if (!exists($counts{$RNA}{$UMI})){
	    # Skip this to save space.
	} else {
	    #	    print STDERR "$RNA\t$UMI\t$counts{$RNA}{$UMI}\n";
	    $sumcount += $counts{$RNA}{$UMI};
	}
    }
    print STD "$RNA\t$sumcount\n";
}

