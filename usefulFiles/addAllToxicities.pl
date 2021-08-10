#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max );

my $toxFile = $ARGV[0];
my $dataFile = $ARGV[1];
my %seeds;
my %viabilities;
my $minLength = 15;

if ($#ARGV < 0){
    die "Usage: addAllToxicities.pl < toxicities file > < data file>\n";
}

my ($mir,$seed,$viability,$header);
open(TOX,$toxFile);
while(<TOX>){
    chomp $_;
    $_ =~ s/\"//g;
    $_ =~ s/\r//g;
    if ($_ !~ /Seed/){
	my($seed,@viabilities) = split(/\t/,$_);
#	$seed =~ s/U/T/g;
	#	print STDERR "$seed\t$viability\n";
	$viabilities{$seed} = "@viabilities";
    } else {
	$header = $_;
    }
}

my (%data);
my $dataHeader = $_;
open(IN,$dataFile);
while(<IN>){
    chomp $_;
    if (($_ !~ /^ID/) && ($_ !~ /^Read/) && ($_ !~ /^\t/)){
	my ($read,@data) = split(/\t/,$_);
	if(length($read) >= $minLength){
	    my $seed = substr($read,1,6);
	    $seed =~ s/T/U/g;
	    if ($seed !~ /N/){
		#	print STDOUT "$id\t$readNum\t$read\t$seed\t$viabilities{$seed}\n";
		my $viabilities = $viabilities{$seed};
		my @viabilities = split(/\s/,$viabilities);		
		my $readData = "$read\t$seed";
		foreach my $viability (@viabilities){
		    $readData .= "\t$viability";
		}
		
		foreach my $datum (@data){
		    $readData .= "\t$datum";
		}
		$data{$read} = $readData;
		$seeds{$read} = $seed;
	    } 
	}
    } else {
	$dataHeader = $_;
    }
}
my @dataHeader = split(/\t/,$dataHeader);
my @viabilityHeader = split(/\t/,$header);

# Print out combined header.
print STDOUT "Read\tSeed";
for (my $i=1;$i<=$#viabilityHeader;$i++){
    print STDOUT "\t$viabilityHeader[$i]";
}
for (my $j=1;$j<=$#dataHeader;$j++){
    print STDOUT "\t$dataHeader[$j]";
}
print STDOUT "\n";
foreach my $read (keys(%data)){
    print STDOUT "$data{$read}\n";
}
