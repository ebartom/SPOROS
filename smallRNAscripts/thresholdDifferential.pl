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

my $prefix = $edgeRfile;
if ($edgeRfile =~ /^([\_\-\.\w\d]+).txt$/){
    $prefix = $1;
    if ($prefix =~ /^([\_\-\.\w\d]+).all$/){
	$prefix = $1;
    }
}

open(UP,">$prefix.up.txt");
#open(DN,">$prefix.dn.txt");
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
	    }
	    if ($labels[$i] eq "logFC"){
		$iLogFC = $i;
		print "Column for LogFC: $iLogFC\n";
	    }
	}
	print UP "$labels[0]";
#	print DN "$labels[0]";
	print DIFF "$labels[0]";
	for (my $i=1;$i<$iLogFC;$i++){
	    print UP "\t$labels[$i]";
#	    print DN "\t$labels[$i]";
	    print DIFF "\t$labels[$i]";
	}
#        print UP "\tNormTotal";
	#        print DN "\tNormTotal";
#	print DIFF "\tNormTotal";
	print UP "\t$labels[-2]";
	print UP "\t$labels[-1]\n";
#	print DN "\t$labels[-2]";
#	print DN "\t$labels[-1]\n";
	print DIFF "\t$labels[-2]";
	print DIFF "\t$labels[-1]\n";
    } else {
	my @data = split(/\t/,$_);
	#	my @countsData = @data[0..$logFC-1];
	my $stat = $data[$iStat];
	if ($stat <= $statThresh){
	    my $logFC = $data[$iLogFC];
#	    print STDERR "$stat\t$statThresh\t$logFC\t$logFCthresh\t$thresholdParameter\n";
	    if ($logFC <= (-1*$logFCthresh)){
#		print DN "$data[0]";
#		for (my $j=1;$j<$iLogFC;$j++){
#		    print DN "\t$data[$j]";
#		}
#		print DN "\t$data[-2]";
#		print DN "\t$data[-1]";
#		print DN "\n";
		print DIFF "$data[0]";
		for (my $j=1;$j<$iLogFC;$j++){
		    print DIFF "\t$data[$j]";
		}
		print DIFF "\t$data[-2]";
		print DIFF "\t$data[-1]";
		print DIFF "\n";
	    } elsif ($logFC >= $logFCthresh){
		print UP "$data[0]";
		for (my $j=1;$j<$iLogFC;$j++){
		    print UP "\t$data[$j]";
		}
		print UP "\t$data[-2]";
		print UP "\t$data[-1]";
		print UP "\n";
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

	    
