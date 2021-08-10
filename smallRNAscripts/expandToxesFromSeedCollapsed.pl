#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::MoreUtils qw(uniq);
use List::Util qw( min max );

my $seedCollapsedFile = $ARGV[0];
my $toxHeader = $ARGV[1];
my $sequenceType = $ARGV[2];
my $project = $ARGV[3];
my $group = $ARGV[4];
my $maxNumSequences = 1000;
my $projectTag = $project;
if ($project ne ""){ $projectTag = "$project.";}

if ($#ARGV <0){
    # This script takes as input the "collapsed" file, a substring
    # to identify the right toxicity / viability column, a sequence
    # type (sRNA / miRNA) and optionally a project ID.
    die "Usage:  $0 <file> <tox type> <sequence type> [project] [group]\n";
}

# Generate some empty variables.
my @toxes;
my %counts;
my @samples;
my %sampleSumCount;
my @headers;
my %toxes;
my $toxIndex = 0;
my $dataIndex = 0;
my @sampleGroup;

my %sampleGroup;
# If a group has been specified, average that, too.
if ($group ne ""){ 
    @sampleGroup = split(/\,/,$group);
    foreach my $sample (@sampleGroup){
	$sampleGroup{$sample} = $project;
    }
    print STDERR "$project\t$group\t@sampleGroup\n";
    $projectTag = "";
}

# Read in the seed collapsed file, one line at a time.
open (IN, $seedCollapsedFile);
while (<IN>){
    chomp $_;
    # If this is the first line, where the first header starts with Seed
    if ($_ =~ /Seed/){
	my $headerLine = $_;
	chomp $headerLine;
	# Then read in the headers.
	@headers = split(/\t/,$headerLine);
	my @headerPrefixes;
	# And iterate through them one by one.
	for (my $i=0;$i<=$#headers;$i++){
	    my $header = $headers[$i];
	    my $lcHeader = lc($header);
	    # When you identify a header that contains the toxHeader
	    # substring, then save its location.
	    # (The tox header can be all lower case and it will still match)
	    if ((($header =~ /$toxHeader/) || ($lcHeader =~ /$toxHeader/)) &&
		($lcHeader =~ /tox/)) {
		$toxIndex = $i;
		print STDERR "$toxHeader Toxicity is in column $toxIndex\n";
	    }
	    # The first five columns are metadata.  After that are data
	    # columns; save the column names.
	    if ($i > 4){
		# Now in data header space.
		$sampleSumCount{$header} = 0;
		if ($group){
		    $sampleSumCount{$project} = 0;
		}
		$header =~ s/\.rep\d//g;
		# If there are multiple replicates of the same sample
		# that end with .rep1 or .rep2, remove the rep suffix
		# and save the sample names as unique header prefixes.
		$sampleSumCount{$header} = 0;
		push(@headerPrefixes,$header);
		#    print "$header\n";
	    }
#	    @samples = uniq (@headerPrefixes);
#	    print STDERR "@samples\n";
	    
	}
	@samples = uniq (@headerPrefixes);
	print STDERR "@samples\n";
    } else {
	# Now done with processing the headers, moving on to the data.
	my @data = split(/\t/,$_);
	# Extract the toxicity / viability data from the right column.
	my $tox = $data[$toxIndex];
	# The miRNA source should be in the second column.
	my $miRNAsource = $data[1];
	# Save the toxicity to a list of all toxicities.
	push (@toxes,$tox);
	if (($sequenceType eq "sRNA") || 
	    ($sequenceType eq "miRNA" && $miRNAsource ne "noMatch")){
	    # If the data is the right sequenceType, then add the counts
	    # to my hash of sample / tox counts.
#	    print STDERR "$sequenceType\t$miRNAsource\n";
	    for (my $i=5;$i<=$#data;$i++){
		if (!exists($counts{$tox}{$headers[$i]})){
		    $counts{$tox}{$headers[$i]} = $data[$i];
		}  else {
		    $counts{$tox}{$headers[$i]} += $data[$i];
		}
		$sampleSumCount{$headers[$i]} += $data[$i];
		if ($group && exists($sampleGroup{$headers[$i]})){
		    $sampleSumCount{$project} += $data[$i];
		}
		# Determine the general name of the sample and add the counts
		# to the general term, too.
		my $sample = $headers[$i];
		$sample =~ s/.rep\d//g;
#		print STDERR "$sample\t$sampleSumCount{$sample}\t$data[$i]\n";
		$sampleSumCount{$sample} += $data[$i];
#		print STDERR "$tox\t$headers[$i]\t$data[$i]\t$counts{$tox}{$headers[$i]}\n";
	    }
	} else {
	    for (my $i=5;$i<=$#data;$i++){
		# If the sequence isn't the right sequence type, still
		# initialize the count for the tox,sample pair.
		if (!exists($counts{$tox}{$headers[$i]})){
		    $counts{$tox}{$headers[$i]} = 0;
		}
	    }
	}
    }
}
# After iterating through every line in the input table, every sample, tox pair
# should be recorded in the counts hash table.



my @replicates;
for (my $i=5;$i<=$#headers;$i++){
    my $header = $headers[$i];
    # Create a list of each data column
    push(@replicates,$header);
    # And print the total counts for each data column.
    print STDERR "$header\t$sampleSumCount{$header}\n";
}
my @sampleSums = values(%sampleSumCount);
my %toxSums;
# Get a unique list of each observed toxicity value.
my @uniqToxes = uniq (@toxes);
#print STDERR "Unique Toxes: @uniqToxes\n";
# Sort the unique toxicities numerically.
@toxes = sort { $a <=> $b }(@uniqToxes);
# Print the unique, sorted toxicities, their number, the samples (without
# replicate suffixes) and the data columns (with replicate suffixes).
# This is all for a sanity check.
print STDERR "Toxes: @toxes\n";
print STDERR "Num Unique Toxes: $#toxes\n";
print STDERR "Samples: @samples\n";
print STDERR "Replicates: @replicates\n";
if ($group){
    print STDERR "Project: $project\n";
    print STDERR "Sample group: @sampleGroup\n";
}
foreach my $sample (@samples){
    print "Sample: $sample\n";
    foreach my $tox (@toxes){
	$toxSums{$tox}{$sample} = 0;
	if ($group && ($group =~ /$sample/)){
	    $toxSums{$tox}{$project} = 0;
	}
    }
    my $numReps = 0;
    foreach my $replicate (@replicates){
#	my $toxRepSum = 0;
	if ($replicate =~ /^$sample/){
	    $numReps++;
	    # For each replicate create an output file.
	    print "Replicate: $replicate\n";
	    my $outfile = "D_toxAnalysis.$replicate.$sequenceType.$toxHeader.".$projectTag."txt";
	    my $normFactor = $sampleSumCount{$replicate}/$maxNumSequences;
#	    print STDERR "$replicate\t$sampleSumCount{$replicate}\t$normFactor\n";
	    open (OUT,">$outfile");
	    foreach my $tox (@toxes){
		# Retrieve the counts for each tox in each data column.
		my $count = $counts{$tox}{$replicate};
		# If replicates of the same sample, sum the counts.
		$toxSums{$tox}{$sample} += $count;
		if ($group && ($group =~ /$sample/)){
		    $toxSums{$tox}{$project} += $count;
		}
#		$toxRepSum+= $count;
		# Print the count as an integer.
		my $intCount = sprintf "%.0f", $count/$normFactor;
		if ($count > 0){
#		    print STDERR "$count\t$intCount\t$normFactor\n";
		}
		if ($intCount > 0){
		    # If the count is greater than 0, print a copy of the
		    # toxicity for each count.
#		    print STDERR "$replicate\t$tox\t$count\t$intCount\t$normFactor\n";
		    for(my $i=1;$i<=$intCount;$i++){
			print OUT "$tox\n";
		    }
		}
		#print STDERR "$toxRepSum\n";
	    }
	    close(OUT);
	    # Print the total counts across the data column.
	    print "ToxCountSum: $sampleSumCount{$replicate}\n";
	}
    }
    if ($numReps > 1){
	# If there is more than one replicate of the data column, also print
	# an average of the replicates.
	my $outfile = "D_toxAnalysis.$sample.avg.$sequenceType.$toxHeader.".$projectTag."txt";
	open(OUT,">$outfile");
	my $normFactor = $sampleSumCount{$sample}/$maxNumSequences;
	print "$sample\t$sampleSumCount{$sample}\tNumReps = $numReps\tnormFactor=$normFactor\n";
	foreach my $tox (@toxes){
	    # Calculate the average count across the replicates and round it.
#	    my $count = ($toxSums{$tox}{$sample})/$numReps;
	    my $count = ($toxSums{$tox}{$sample})/$normFactor;
	    my $intCount = sprintf "%.0f", $count;
#	    if ($sample eq "DicerKO"){
#		print STDERR "$sample\t$tox\t$toxSums{$tox}{$sample}\t$numReps\t$count\t$intCount\n";
#	    }
	    if ($intCount > 0){
		# For each count, print a copy of the toxicity.
		for (my $i=1;$i<=$intCount;$i++){
		    print OUT "$tox\n";
		}
	    }
	}
	close(OUT);
    }
}
if ($group){
    # If there is a defined sample group, also print
    # an average of the replicates.
    my $outfile = "D_toxAnalysis.$project.avg.$sequenceType.$toxHeader.".$projectTag."txt";
    open(OUT,">$outfile");
    my $normFactor = $sampleSumCount{$project}/$maxNumSequences;
    print "$project\t$sampleSumCount{$project}\tNumReps = ".($#sampleGroup+1)."\tnormFactor=$normFactor\n";
    foreach my $tox (@toxes){
	    # Calculate the average count across the replicates and round it.
	#	    my $count = ($toxSums{$tox}{$sample})/$numReps;
	my $count = ($toxSums{$tox}{$project})/$normFactor;
	    my $intCount = sprintf "%.0f", $count;
	#	    if ($sample eq "DicerKO"){
	#		print STDERR "$sample\t$tox\t$toxSums{$tox}{$sample}\t$numReps\t$count\t$intCount\n";
	#	    }
	if ($intCount > 0){
	    # For each count, print a copy of the toxicity.
	    for (my $i=1;$i<=$intCount;$i++){
		print OUT "$tox\n";
	    }
	}
    }
    close(OUT);
}

