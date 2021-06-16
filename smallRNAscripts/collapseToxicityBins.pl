#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $collapsedSpeciesFile = $ARGV[0];
my $toxHeader = $ARGV[1];
my $sequenceType = $ARGV[2];

if ($#ARGV <0){
    die "Usage:  $0 <file> <tox type> <sequence type>\n";
}

if ($sequenceType eq ""){
    $sequenceType = "sRNA";
}
print STDERR "SequenceType = $sequenceType\n";

my $prefix = basename($collapsedSpeciesFile);
if ($prefix =~ /^([.\.\-\_\w]+).txt$/){
    $prefix = $1;
    # Remove the "collapsed" part to shorten the filenames.
    $prefix =~ s/collapsed.//g;
}
print STDERR "$prefix\n";
my $outfile1 = "seedKeyed.$sequenceType.$toxHeader.$prefix.txt";
my $outfile2 = "binned.$sequenceType.$toxHeader.$prefix.txt";

open (IN,$collapsedSpeciesFile);;
open (OUT1,">$outfile1");
open (OUT2,">$outfile2");

my %seeds;
my @headers;
my $toxIndex = 0;
my $dataIndex = 0;
my $sequenceTypeFlag = 0;
my %toxes;
my %toxesKeyed;
my ($seed, $tox, @data);
my ($source1, $source2, $tox1, $tox2);
my $toxLabel = "";
my $headerLine = "";
print STDERR "Tox data = $toxHeader\n";

while(<IN>){
    chomp $_;
    # First find the right tox column (even if lower case or upper case)
    if ($_ =~ /^Seed/){
	$headerLine = $_;
	@headers = split(/\t/,$headerLine);
	$dataIndex = 0;
	my $lcToxHeader = lc $toxHeader;
	foreach my $header (@headers){
	    #	    print STDERR "$header\n";
	    my $lcHeader = lc $header;
	    if ($lcHeader =~ /$lcToxHeader/){
#		print "$toxHeader\t$header\n";
		$toxIndex = $dataIndex;
		$toxLabel = $header;
	    }
	    $dataIndex++;
	}
    }else{
	$sequenceTypeFlag = 0;
	($seed,@data) = split(/\t/,$_);
	$tox = $data[$toxIndex-1];
	($seed,$source1,$source2,$tox1,$tox2,@data) = split(/\t/,$_);
	if ($sequenceType eq "sRNA"){
	    # Include Everything.
	    $sequenceTypeFlag = 1;
	} elsif ($sequenceType eq "miRNA"){
	    # Only include if there is a hit in source1 (miRNA blast result).
	    if ($source1 ne "noMatch"){
		$sequenceTypeFlag = 1;
	    }
	}
	if ($sequenceTypeFlag == 1){
	    # Include this row in analysis.
	    if (!exists($seeds{$seed})){
		$seeds{$seed} = "@data";
		$toxes{$seed} = $tox;
	    } else {
		my @countsOld = split(/\s/,$seeds{$seed});
		my @countsNew = @data;
		for (my $i=0;$i<=$#countsNew;$i++){
		    $countsNew[$i] +=$countsOld[$i];
		}
		$seeds{$seed} = "@countsNew";
	    }
	    my $roundedTox = sprintf "%.0f", $tox;
	    if (!exists($toxesKeyed{$roundedTox})){
		$toxesKeyed{$roundedTox} = "@data";
	    }else {
	    my @countsOld = split(/\s/,$toxesKeyed{$roundedTox});
	    my @countsNew = @data;
	    for (my $i=0;$i<=$#countsNew;$i++){
		$countsNew[$i] +=$countsOld[$i];
	    }
	    $toxesKeyed{$roundedTox} = "@countsNew";
	    }
	}
    }
}
close(IN);

my @sortedSeeds = sort(keys(%seeds));
print OUT1 "Seed\t$toxLabel";
my ($seedHeader,$source1header,$source2header,$tox1header,$tox2header,@dataHeaders) = split(/\t/,$headerLine);
#print STDERR @dataHeaders;
foreach my $dataHeader (@dataHeaders){
    print OUT1 "\t$dataHeader";
}
print OUT1 "\n";

foreach my $seed (@sortedSeeds){
    print OUT1 "$seed\t$toxes{$seed}";
    my @counts = split(/\s/,$seeds{$seed});
    foreach my $count (@counts){
	print OUT1 "\t$count";
    }
    print OUT1 "\n";
}
close(OUT1);

print OUT2 $toxLabel;
foreach my $dataHeader (@dataHeaders){
    print OUT2 "\t$dataHeader";
}
print OUT2 "\n";
my @sortedToxes = sort {$a <=> $b} (keys(%toxesKeyed));
foreach my $tox (@sortedToxes){
    print OUT2 "$tox";
    my @counts = split(/\s/,$toxesKeyed{$tox});
    foreach my $count (@counts){
	print OUT2 "\t$count";
    }
    print OUT2 "\n";
}
close(OUT2);
    

