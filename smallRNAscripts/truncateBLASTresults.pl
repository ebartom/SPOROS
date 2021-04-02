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

while(<AN>){
    chomp $_;
    my ($read,@data) = split(/\t/,$_);
    print ($read);
    foreach my $datum (@data){
	if ( $datum =~ /\:/){
#	    print STDERR "Before: $datum\n";
	    if ($datum =~ /^([\.\-\_\|\w]+)\:/){
		$datum = $1;
		print "\t$datum";
#		print STDERR "After: $datum\n";
	    }
	} else {
#	    print STDERR "\nNOT: $datum\n";
	    print "\t$datum";
	}
    }
    print "\n";
}
