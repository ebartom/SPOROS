#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $originalfile = $ARGV[0];
my $blastfile = $ARGV[1];

my %blastresults;

if ($#ARGV <0){
    die "Usage:  $0 <file> <blastfile>\n";
}


open (BL,$blastfile);
print STDERR "Blast file #1: $blastfile\n";

# Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score

while(<BL>){
    chomp $_;
    if ($_ !~ /^#/){
	my($seq,$subject,$perID,$length,$mism,$gapOpens,$qstart,$qend,$sstart,$send,$evalue,$bitscore)= split(/\t/,$_);
#	print STDERR "$subject\n";
	$seq =~ s/\.\d+//;
#	print STDERR "BL \"$seq\"\n";
	if (!exists($blastresults{$seq})){
	    $blastresults{$seq} = "$subject\:$qstart-$qend\:$sstart-$send\:$perID\:$evalue";
	} else {
	    $blastresults{$seq} .= "|$subject\:$qstart-$qend\:$sstart-$send\:$perID\:$evalue";
	}
    }
}
close(BL);


open (IN,$originalfile);
while(<IN>){
    chomp $_;
    if (($_ !~ /^Read/) && ($_ !~ /^\t/)){
	my ($seq,@data) = split(/\t/,$_);
#    print STDERR "IN: \"$seq\"\n";
#    my $id = $data[-1];
#    print STDERR "ID:  $id\n";
#    my $label = "$seq.$id";
	if (exists ($blastresults{$seq})){
	    print "$_\t$blastresults{$seq}\n";
	} else {
	    print "$_\tnoMatch\n";
	}
    } else {
	print "$_\t$blastfile\n";
    }
}
