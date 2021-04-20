#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $annotatedfile = $ARGV[0];

if ($#ARGV <0){
    die "Usage:  $0 <file>\n";
}

open (AN,$annotatedfile);

my $numSamples = 0;
my @dataIndices;
my @headers;
my %allSpecies;
my %toxes;
my $dataIndex = 0;
my $toxIndex = 0;
my $labelIndex = 0;
my($read,$seed,@data,$label);
while(<AN>){
    chomp $_;
    if ($_ =~ /^Read/){
	@headers = split(/\t/,$_);
	$dataIndex = 0;
	foreach my $header (@headers){
	    if ($header =~ /uniqCounts/){
		push(@dataIndices,$dataIndex);
		$numSamples++;
	    }
	    if ($header =~ /RNAworld/){
		$labelIndex = $dataIndex;
	    }
	    if ($header =~ /3human/){
		$toxIndex = $dataIndex;
	    }
	    $dataIndex++;
	}
#	print STDERR "@dataIndices\n";
    }else{
	($read,$seed,@data) = split(/\t/,$_);
	#print ($read);
       $label = $data[$labelIndex-2];
	if ($label =~ /^([\.\-\_\w]+)\|?/){
	    $label = $1;
	}
	my $species = "$seed\t$label";
	my $tox = $data[$toxIndex-2];
	$toxes{$species} = $tox;
	my @counts;
	foreach my $index (@dataIndices){
	    push (@counts,$data[$index-2]);
	}
	if (!exists($allSpecies{$species})){
	    $allSpecies{$species} = "@counts";
	} else {
	    my @counts2 = split(/\s/,$allSpecies{$species});
	    #	print STDERR "$species\t@counts\t@counts2\	
	    for (my $index=0;$index<=$#counts;$index++){
		$counts2[$index] += $counts[$index];
	    }
	    #	print STDERR "$species\t@counts2\n";
	    $allSpecies{$species} = "@counts2";
	}
    }
}

print "Seed\tRNAworld BlastHit";
my $toxLabel = $headers[$toxIndex];
print "\t$toxLabel";

foreach my $dataIndex (@dataIndices){
    my $header = $headers[$dataIndex];
    $header =~ s/.justReads.uniqCounts.txt//g;
    print "\t$header";
}

print "\n";
foreach my $seq (keys(%allSpecies)){
    print "$seq";
    my $tox = $toxes{$seq};
    print "\t$tox";
    my @counts = split(/\s/,$allSpecies{$seq});
    foreach my $count (@counts){
	print "\t$count";
    }
    print "\n";
}
