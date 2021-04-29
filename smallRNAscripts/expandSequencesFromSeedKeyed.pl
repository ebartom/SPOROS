#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

my $seedKeyedFile = $ARGV[0];
my $maxNumSequences = 1000;

if ($#ARGV <0){
    die "Usage:  $0 <file>\n";
}

my @seeds;
my %counts;
my @samples;
my %sampleSumCount;
my @headers;
open (IN, $seedKeyedFile);
while (<IN>){
    chomp $_;
    if ($_ =~ /^Seed/){
	my $headerLine = $_;
	chomp $headerLine;
	@headers = split(/\t/,$headerLine);
	my @headerPrefixes;
	for (my $i=2;$i<=$#headers;$i++){
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
	my $seed = $data[0];
	push (@seeds,$seed);
	for (my $i=2;$i<=$#data;$i++){
	    $counts{$seed}{$headers[$i]} = $data[$i];
	    $sampleSumCount{$headers[$i]} += $data[$i];
#	    print "$seed\t$headers[$i]\t$data[$i]\n";
	}
    }
}

my @replicates;
for (my $i=2;$i<=$#headers;$i++){
    my $header = $headers[$i];
    push(@replicates,$header);
    print "$header\t$sampleSumCount{$header}\n";
}
my @sampleSums = values(%sampleSumCount);
my $maxCounts = max(@sampleSums);
print "$maxCounts\n";
my $normFactor = $maxCounts/$maxNumSequences;
print "$normFactor\n";

my %seedSums;
foreach my $sample (@samples){
    print "Sample: $sample\n";
    foreach my $seed (@seeds){
	$seedSums{$seed}{$sample} = 0;
    }
    my $numReps = 0;
    foreach my $replicate (@replicates){
	if ($replicate =~ /^$sample/){
	    $numReps++;
	    print "Replicate: $replicate\n";
	    my $outfile = "$replicate.seeds.$maxNumSequences.txt";
	    open (OUT,">$outfile");
#	    print OUT ">$replicate.seeds $maxNumSequences\n";
	    my $seedRepSum = 0;
	    foreach my $seed (@seeds){
		my $count = $counts{$seed}{$replicate};
		my $normCount = sprintf "%.0f", ($count/$normFactor);
		$seedSums{$seed}{$sample} += $normCount;
		$seedRepSum+= $normCount;
		if ($normCount > 0){
#		    print "$seed\t$normCount\n";
		    for(my $i=0;$i<=$normCount;$i++){
			print OUT "$seed\n";
		    }
		}
	    }
	    close(OUT);
	    print "SeedSum: $seedRepSum\n";
	}
    }
    my $outfile = "$sample.avg.seeds.$maxNumSequences.txt";
    open(OUT,">$outfile");
#    print OUT ">$sample.avg.seeds $maxNumSequences\n";
    print "$sample\tNumReps = $numReps\n";
    foreach my $seed (@seeds){
	my $count = $seedSums{$seed}{$sample}/$numReps;
	my $intCount = sprintf "%.0f", $count;
	if ($intCount > 0){
	    for (my $i=0;$i<=$intCount;$i++){
		print OUT "$seed\n";
	    }
	}
    }
    close(OUT);
}
