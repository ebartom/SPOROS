#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $collapsedSpeciesFile = $ARGV[0];
my $toxHeader = $ARGV[1];

if ($#ARGV <0){
    die "Usage:  $0 <file> <tox type>\n";
}

open (IN,$collapsedSpeciesFile);;

my %seeds;
my @headers;
my $toxIndex = 0;
my $dataIndex = 0;
my %toxes;
my ($seed, $tox, @data);
my ($source1, $source2, $tox1, $tox2);
my $toxLabel = "";
my $headerLine = "";
print STDERR "Tox data = $toxHeader\n";

while(<IN>){
    chomp $_;
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
	($seed,@data) = split(/\t/,$_);
	$tox = $data[$toxIndex-1];
	($seed,$source1,$source2,$tox1,$tox2,@data) = split(/\t/,$_);
	if (!exists($seeds{$seed})){
	    $seeds{$seed} = "@data";
	    $toxes{$seed} = $tox;
	} else {
	    my @countsOld = split(/\s/,$seeds{$seed});
	    for (my $i=0;$i<=$#data;$i++){
		$data[$i] +=$countsOld[$i];
	    }
	    $seeds{$seed} = "@data";
	}
    }
}

my @sortedSeeds = sort(keys(%seeds));
print "Seed\t$toxLabel";
my ($seedHeader,$source1header,$source2header,$tox1header,$tox2header,@dataHeaders) = split(/\t/,$headerLine);
foreach my $dataHeader (@dataHeaders){
    print "\t$dataHeader";
}
print "\n";

foreach my $seed (@sortedSeeds){
    print "$seed\t$toxes{$seed}";
    my @counts = split(/\s/,$seeds{$seed});
    foreach my $count (@counts){
	print "\t$count";
    }
    print "\n";
}
