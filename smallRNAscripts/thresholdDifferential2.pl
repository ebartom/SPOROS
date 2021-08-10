#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $edgeRfile          = $ARGV[0];
my $thresholdParameter = $ARGV[1];
my $logFCthresh        = $ARGV[2];
my $statThresh         = $ARGV[3];

if ($#ARGV <0){
    die "Usage:  $0 <file> <parameter> <logFCthresth> <statThresh>\n";
}

if ($thresholdParameter eq "adjp"){
    $thresholdParameter = "adj.p";
}
if (($thresholdParameter eq "PValue") || ($thresholdParameter eq "pvalue")){
    $thresholdParameter = "PValue";
}

my $iStat = 0;
my $iLogFC = 0;

my $prefix = basename($edgeRfile);
if ($prefix =~ /^([\_\-\.\w\d]+).txt$/){
    $prefix = $1;
    if ($prefix =~ /^([\_\-\.\w\d]+).all$/){
	$prefix = $1;
    }
}

open(DIFF,">$prefix.diff.txt");

open(IN,$edgeRfile);
while(<IN>){
    chomp $_;
    if (($_ =~ /^Read/) ){
	print "$_\n";
	my @labels = split(/\t/,$_);
	for(my $i=0;$i<=$#labels;$i++){
	    if($labels[$i] eq $thresholdParameter){
		$iStat = $i;
		print "Column for $thresholdParameter: $iStat\n";
	    }
	    if ($labels[$i] eq "logFC"){
		$iLogFC = $i;
		print "Column for LogFC: $iLogFC\n";
	    }
	}
	print DIFF "$labels[0]";
	for (my $i=1;$i<$iLogFC;$i++){
	    print DIFF "\t$labels[$i]";
	}
	print DIFF "\t$labels[-2]";
	print DIFF "\t$labels[-1]\n";
    } else {
	my @data = split(/\t/,$_);
	my $stat = $data[$iStat];
	if ($stat <= $statThresh){
	    my $logFC = $data[$iLogFC];
#	    print STDERR "$stat\t$statThresh\t$logFC\t$logFCthresh\t$thresholdParameter\n";
	    if ($logFC <= (-1*$logFCthresh)){
		print DIFF "$data[0]";
		for (my $j=1;$j<$iLogFC;$j++){
		    print DIFF "\t$data[$j]";
		}
		print DIFF "\t$data[-2]";
		print DIFF "\t$data[-1]";
		print DIFF "\n";
	    } elsif ($logFC >= $logFCthresh){
		print DIFF "$data[0]";
		for (my $j=1;$j<$iLogFC;$j++){
		    print DIFF "\t$data[$j]";
		}
		print DIFF "\t$data[-2]";
		print DIFF "\t$data[-1]";
		print DIFF "\n";
	    }
	}
    }
}

	    
