#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $originalfile = $ARGV[0];

if ($#ARGV <0){
    die "Usage:  $0 <file>\n";
}



open (IN,$originalfile);
while(<IN>){
    chomp $_;
    if (($_ !~ /^Read/) && ($_ !~ /^\t/)){
	my ($seq,@data) = split(/\t/,$_);
#	print STDERR "$seq\n";
	if ($seq =~ /GTCCGACGATC/){
#	    	if ($seq =~ /GTCCGACGATC[ACTG]{3,5}$/){
#	    print STDERR "Found adapter: $seq\n";
	} else {
	    # This sequence does not look like adapter.
	    print "$_\n";
	}
    } else {
	# This is the header line.
	print "$_\n";
    }
}
