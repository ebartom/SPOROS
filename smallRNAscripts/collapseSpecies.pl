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
my $humanToxIndex = 0;
my $mouseToxIndex = 0;
my $label_miRNAsIndex = 0;
my $label_RNAworldIndex = 0;
my($read,$seed,@data,$miRNAlabel,$RNAworldlabel);
while(<AN>){
    chomp $_;
    if ($_ =~ /^Read/){
	@headers = split(/\t/,$_);
	$dataIndex = 0;
	foreach my $header (@headers){
	    if (($header =~ /uniqCounts/) || ($header=~ /rep\d/) || ($header =~ /SRR/)){
		push(@dataIndices,$dataIndex);
		$numSamples++;
	    }
	    if ($header =~ /miRNAs.nr/){
		$label_miRNAsIndex = $dataIndex;
	    }
	    if ($header =~ /RNAworld/){
		$label_RNAworldIndex = $dataIndex;
	    }
	    if ($header =~ /3human/){
		$humanToxIndex = $dataIndex;
	    }
	    if ($header =~ /3mouse/){
		$mouseToxIndex = $dataIndex;
	    }
	    $dataIndex++;
	}
#	print STDERR "@dataIndices\n";
    }else{
	($read,$seed,@data) = split(/\t/,$_);
	#print ($read);
       $miRNAlabel = $data[$label_miRNAsIndex-2];
	if ($miRNAlabel =~ /^([\.\-\_\w]+)\|?/){
	    $miRNAlabel = $1;
	}
	$RNAworldlabel = $data[$label_RNAworldIndex-2];
	if ($RNAworldlabel =~ /^([\.\-\_\w]+)\|?/){
	    $RNAworldlabel = $1;
	}
	my $species = "$seed\t$miRNAlabel\t$RNAworldlabel";
	my $toxHuman = $data[$humanToxIndex-2];
	my $toxMouse = $data[$mouseToxIndex-2];
	$toxes{$species} = "$toxHuman\t$toxMouse";
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

print "Seed\tTop miRNA Source\tTop RNAworld Source";
#my $toxLabel = $headers[$toxIndex];
print "\tAvg of 3 Human Tox\tAvg of 3 Mouse Tox";

foreach my $dataIndex (@dataIndices){
    my $header = $headers[$dataIndex];
    $header =~ s/.justReads.uniqCounts.txt//g;
    $header =~ s/.justReads.uniqCounts.notUMId.txt//g;
    $header =~ s/.justReads.uniqCounts.UMId.txt//g;
    $header =~ s/.noAdapter.txt//g;
    print "\t$header";
}

print "\n";
my @sortedRNAspecies = sort (keys(%allSpecies));
foreach my $seq (@sortedRNAspecies){
    print "$seq";
    my $tox = $toxes{$seq};
    print "\t$tox";
    my @counts = split(/\s/,$allSpecies{$seq});
    foreach my $count (@counts){
	print "\t$count";
    }
    print "\n";
}
