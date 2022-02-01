#!/usr/bin/perl -w
use strict;
use warnings;

my $rawCountsTable = $ARGV[0];
my $minSum = $ARGV[1];
if ($minSum) {
    # MinSum was defined on the command line.
}else {
    $minSum = 2;
}
print STDERR "Remove reads with fewer than $minSum counts across all samples.\n";

my %counts;
my %countTotals;
my %reads;
my @headers;
my $numSamples = 0;
my $stuff = "";

open(IN,$rawCountsTable);
while (<IN>){
    chomp $_;
    if (($_ =~ /^Read/) || ($_ =~ /^Symbol/)){
	($stuff,@headers) = split(/\t/,$_);
	$numSamples = $#headers;
	if (($headers[$numSamples] eq "NormTotal") || 
	    ($headers[$numSamples] eq "RowTotal")){
	    $numSamples--;
	}
	for (my $i=0; $i<= $numSamples;$i++){
	    $countTotals{$i}= 0;
	}
    } else {
	my ($read,@data) = split(/\t/,$_);
	$reads{$read} = 0;
	for (my $i=0; $i <= $numSamples;$i++){
	    my $count = $data[$i];
	    $counts{$read}{$i} = $count;
	    $countTotals{$i}+= $count;
	    $reads{$read} += $count;
	}
    }
}

print STDERR "Total number of unique reads\n";
my @reads = keys(%reads);
print STDERR scalar(@reads);
print "Read";
for (my $i=0;$i<=$numSamples;$i++){
    print "\t$headers[$i]";
}
print "\tNormTotal\n";
foreach my $read (keys(%reads)){
    if ($reads{$read} >= $minSum){
	print $read;
	my $total = 0;
	for (my $i=0;$i<=$numSamples;$i++){
	    if (exists($counts{$read}{$i})){
		my $normCount = ($counts{$read}{$i}/$countTotals{$i})*1000000;
		printf "\t%.3f",$normCount;
		$total += $normCount;
#		print "\t$counts{$read}{$i}";
	    } else {
		print "\t0";
	    }
	}
	print "\t$total";
	print "\n";
    }
}
