#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

my $seedKeyedFile = $ARGV[0];
my $sequenceType = $ARGV[1];
my $project = $ARGV[2];
my $maxNumSequences = 1000;
my $seedLength = 6;

if ($project ne ""){ $project = "$project.";}

if ($seedKeyedFile !~ /$sequenceType/){
    print "WARNING: This input file does not seem specific to $sequenceType!\n";
}
    

if ($#ARGV <0){
    die "Usage:  $0 <file> <sequence type> [project]\n";
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
#my $maxCounts = max(@sampleSums);
#print "$maxCounts\n";
#my $normFactor = $maxCounts/$maxNumSequences;
#print "$normFactor\n";

my %seedSums;
foreach my $sample (@samples){
    print "Sample: $sample\n";
    foreach my $seed (@seeds){
	$seedSums{$seed}{$sample} = 0;
    }
    my $numReps = 0;
    for (my $i=0;$i<=$#replicates;$i++){
	my $replicate = $replicates[$i];
#    foreach my $replicate (@replicates){
	if ($replicate =~ /^$sample/){
	    $numReps++;
	    print "Replicate: $replicate\n";
	    my $outfile1 = "seedAnalysis.$replicate.$sequenceType.".$project."$maxNumSequences.txt";	    
	    my $outfile2 = "seedLong.$replicate.$sequenceType.".$project."$maxNumSequences.txt";
	    open (OUT1,">$outfile1");
	    open (OUT2,">$outfile2");
	    #	    print OUT ">$replicate.seeds $maxNumSequences\n";
	    print OUT2 "Seed\tSample\tSeedID\tPos\tBase\n";
	    my $seedRepSum = 0;
	    my $seedNum = 0;
	    my $sampleSumCount = $sampleSumCount{$replicate};
	    if (exists($sampleSumCount{$sample})){
		$sampleSumCount{$sample}+=$sampleSumCount{$replicate};
	    } else {
		$sampleSumCount{$sample} = $sampleSumCount{$replicate};
	    }
	    my $normFactor = $sampleSumCount/$maxNumSequences;
	    print STDERR "$replicate \t$normFactor\n";
	    foreach my $seed (@seeds){
		my $count = $counts{$seed}{$replicate};
		my $normCount = sprintf "%.0f", ($count/$normFactor);
		$seedSums{$seed}{$sample} += $normCount;
		$seedRepSum+= $normCount;
		if ($normCount > 0){
#		    print "$seed\t$normCount\n";
		    for(my $i=1;$i<=$normCount;$i++){
			$seedNum++;
			print OUT1 "$seed\n";
			for (my $j=1;$j<=$seedLength;$j++){
			    my $base = substr($seed,$j-1,1);
			    print OUT2 "$seed\t$replicate\t$replicate.$seedNum\t$j\t$base\n";
			}
		    }
		}
	    }
	    close(OUT1);
	    close(OUT2);
	    print "SeedSum: $seedRepSum\n";
	}
    }
    if ($numReps > 1){
	my $outfile1 = "seedAnalysis.$sample.avg.$sequenceType.".$project."$maxNumSequences.txt";
	my $outfile2 = "seedLong.$sample.avg.$sequenceType.".$project."$maxNumSequences.txt";
	open(OUT1,">$outfile1");
	open(OUT2,">$outfile2");
	print OUT2 "Seed\tSample\tSeedID\tPos\tBase\n";
	print STDERR "$sample\tNumReps = $numReps\n";
	my $seedNum=0;
	foreach my $seed (@seeds){
	    my $count = $seedSums{$seed}{$sample}/$numReps;
	    my $intCount = sprintf "%.0f", $count;
	    if ($intCount > 0){
		for (my $i=1;$i<=$intCount;$i++){
		    $seedNum++;
		    print OUT1 "$seed\n";
		    for (my $j=1;$j<=$seedLength;$j++){
			my $base = substr($seed,$j-1,1);
			print OUT2 "$seed\t$sample.avg\t$sample.$seedNum\t$j\t$base\n";
		    }
		}
	    }
	}
	close(OUT1);
	close(OUT2);
    }
}
