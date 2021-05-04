#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

my $toxKeyedFile = $ARGV[0];
my $maxNumToxes = 1000;

if ($#ARGV <0){
    die "Usage:  $0 <file>\n";
}

my @toxes;
my %counts;
my @samples;
my %sampleSumCount;
my @headers;
open (IN, $toxKeyedFile);
while (<IN>){
    chomp $_;
    if ($_ =~ /Tox/){
	my $headerLine = $_;
	chomp $headerLine;
	@headers = split(/\t/,$headerLine);
	my @headerPrefixes;
	for (my $i=1;$i<=$#headers;$i++){
	    my $header = $headers[$i];
	    $sampleSumCount{$header} = 0;
	    $header =~ s/\.rep\d//g;
	    push(@headerPrefixes,$header);
	    #    print "$header\n";
	}
	@samples = uniq (@headerPrefixes);
	print "@samples\n";
    } else {
	my @data = split(/\t/,$_);
	my $tox = $data[0];
	push (@toxes,$tox);
	for (my $i=1;$i<=$#data;$i++){
	    $counts{$tox}{$headers[$i]} = $data[$i];
	    $sampleSumCount{$headers[$i]} += $data[$i];
#	    print "$tox\t$headers[$i]\t$data[$i]\n";
	}
    }
}

my @replicates;
for (my $i=1;$i<=$#headers;$i++){
    my $header = $headers[$i];
    push(@replicates,$header);
    print "$header\t$sampleSumCount{$header}\n";
}
my @sampleSums = values(%sampleSumCount);
my $maxCounts = max(@sampleSums);
print "$maxCounts\n";
my $normFactor = $maxCounts/$maxNumToxes;
print "$normFactor\n";

my %toxSums;
foreach my $sample (@samples){
    print "Sample: $sample\n";
    foreach my $tox (@toxes){
	$toxSums{$tox}{$sample} = 0;
    }
    my $numReps = 0;
    foreach my $replicate (@replicates){
	if ($replicate =~ /^$sample/){
	    $numReps++;
	    print "Replicate: $replicate\n";
	    my $outfile = "$replicate.toxes.$maxNumToxes.txt";
	    open (OUT,">$outfile");
#	    print OUT ">$replicate.toxes $maxNumToxes\n";
	    my $toxRepSum = 0;
	    foreach my $tox (@toxes){
		my $count = $counts{$tox}{$replicate};
		my $normCount = sprintf "%.0f", ($count/$normFactor);
		$toxSums{$tox}{$sample} += $normCount;
		$toxRepSum+= $normCount;
		if ($normCount > 0){
#		    print "$tox\t$normCount\n";
		    for(my $i=0;$i<=$normCount;$i++){
			print OUT "$tox\n";
		    }
		}
	    }
	    close(OUT);
	    print "ToxSum: $toxRepSum\n";
	}
    }
    my $outfile = "$sample.avg.toxes.$maxNumToxes.txt";
    open(OUT,">$outfile");
    print "$sample\tNumReps = $numReps\n";
    foreach my $tox (@toxes){
	my $count = $toxSums{$tox}{$sample}/$numReps;
	my $intCount = sprintf "%.0f", $count;
	if ($intCount > 0){
	    for (my $i=0;$i<=$intCount;$i++){
		print OUT "$tox\n";
	    }
	}
    }
    close(OUT);
}
