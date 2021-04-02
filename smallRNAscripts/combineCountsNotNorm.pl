#!/usr/bin/perl -w
use strict;
use warnings;

my $fileSuffix = $ARGV[0];
if ($fileSuffix) {
    print STDERR "$fileSuffix\n";
    # Great!
} else {
    $fileSuffix = "justReads.uniqCounts.txt";
    print STDERR "Default suffix: $fileSuffix\n";
}
print STDERR "Look for files with Suffix $fileSuffix\n";
my $countfiles;
$countfiles = `ls *.$fileSuffix`;
my @countfiles = split(/\s+/,$countfiles);
my $numFiles = scalar(@countfiles);
my $minSum = 2;
if ($fileSuffix =~ /\.UMId/){$minSum = 1;}
print STDERR "MinSum: $minSum\n";
print STDERR "Countfiles:\n@countfiles\n";

#my $umiCheck = $ARGV[0];
#my $countfiles;
#print STDERR "Using UMId files? $umiCheck\n";
#if ($umiCheck == 1){
#    $countfiles = `ls *.uniqCounts.UMId.txt`;
#}else{
#    $countfiles = `ls *.uniqCounts.txt`;
#}
#my @countfiles = split(/\s+/,$countfiles);
#my $numFiles = scalar(@countfiles);
#my $minSum = 2;
#if ($umiCheck ==1) {$minSum = 1;};
#print STDERR "MinSum: $minSum\n";
#print STDERR "Countfiles:\n@countfiles\n";

my %counts;
my %countTotals;
my %reads;
my $fileindex = 0;
foreach my $file (@countfiles){
    print STDERR "$file\n";
    $fileindex++;
    $countTotals{$fileindex}= 0;
    open(IN,$file);
    while (<IN>){
	my ($stuff,$count,$read) = "";
	if ($fileSuffix =~ /UMId/){
#	if ($umiCheck == 1){
	    ($read,$count) = split(/\s+/,$_);
	} else {
	    ($stuff,$count,$read) = split(/\s+/,$_);
	}
	if ($stuff ne ""){ $read = $count; $count = $stuff;}
#	print "$read,$count\n";
	if (exists ($reads{$read})){
	    $reads{$read} += $count;
	} else {$reads{$read} = $count;}
	$counts{$read}{$fileindex} = $count;
	$countTotals{$fileindex}+= $count;
    }
}

print STDERR "Total number of unique reads\n";
my @reads = keys(%reads);
print STDERR scalar(@reads);
print "Read";
foreach my $file (@countfiles){
    print "\t$file";
}
print "\tNormTotal\n";
foreach my $read (keys(%reads)){
    if ($reads{$read} >= $minSum){
	print $read;
	my $total = 0;
	for (my $i=1;$i<=$fileindex;$i++){
	    if (exists($counts{$read}{$i})){
		#my $normCount = ($counts{$read}{$i}/$countTotals{$i})*1000000;
#		printf "\t%.3f",$normCount;
		$total += $counts{$read}{$i};
		print "\t$counts{$read}{$i}";
	    } else {
		print "\t0";
	    }
	}
	print "\t$total";
	print "\n";
    }
}
