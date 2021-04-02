#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $originalfile = $ARGV[0];
my $threshold = $ARGV[1];


if ($#ARGV <0){
    die "Usage:  $0 <file> <threshold>\n";
}

my $iSum = 0;
open(IN,$originalfile);
while(<IN>){
    chomp $_;
    if (($_ =~ /^Read/) ){
	print "$_\n";
	my @labels = split(/\t/,$_);
	for(my $i=0;$i<=$#labels;$i++){
	    if($labels[$i] eq "NormTotal"){
		$iSum = $i;
	    }
	}
    } else {
	my @data = split(/\t/,$_);
	if ($data[$iSum] >= $threshold){
	    print "$_\n";
	}
    }
}

	    
