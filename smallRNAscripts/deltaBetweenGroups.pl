#!/opt/local/bin/perl -w

use strict;
use File::Basename;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < dataTable >\n";
  exit;
}

my $dataTable = $ARGV[0];
my $prefix = basename($dataTable);
if ($prefix =~ /^([\_\-\.\w\d]+).txt$/){
    $prefix = $1;
    if ($prefix =~ /^([\_\-\.\w\d]+).diff.avgd$/){
	$prefix = $1;
    }
}
open(UP,">$prefix.delta.up.txt");
open(DN,">$prefix.delta.dn.txt");


my $group1Col = 0;
my $group2Col = 0;
my (@labels);
open (IN,$dataTable);
while (<IN>){
    chomp $_;
    if (($_ =~ /^Read/) || ($_ =~ /^Avg/) || ($_ =~ /^Seed/)){
	@labels = split(/\t/,$_);
	for (my $i=0;$i<=$#labels;$i++){
	    my $label = $labels[$i];
	    if ($label =~ /Group1/){
		$group1Col = $i;
	    }
	    if ($label =~ /Group2/){
		$group2Col = $i;
	    }
	}
	print UP "$labels[0]\t$labels[1]\tDelta\n";
	print DN "$labels[0]\t$labels[1]\tDelta\n";
    } else {
	my @data = split(/\t/,$_);
	my $delta = $data[$group2Col] - $data[$group1Col];
	if ($delta > 0){
	    print UP "$data[0]\t$data[1]\t$delta\n";
	} elsif ($delta < 0){
	    print DN "$data[0]\t$data[1]\t".-1*$delta."\n";
	}
    }
}
