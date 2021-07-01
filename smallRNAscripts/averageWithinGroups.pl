#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < comparisonsFile > < dataTable > < comparison >\n";
  exit;
}

my $compfile = $ARGV[0];
my $dataTable = $ARGV[1];
my $myComparison = $ARGV[2];

my @labels;
my @data;
my %samples;
my %group1samples;
my %group2samples;

open (CMP,$compfile);
while (<CMP>){
    chomp $_;
    if ($_ =~ /^Comparison/){
	@labels = split(/\,/,$_);
	for (my $i=0;$i<=$#labels;$i++){
	    $samples{$i} = $labels[$i];
#	    print "$labels[$i]\n";
	}
    } else {
	@data = split(/\,/,$_);
	my $comparison = $data[0];
#	print STDERR "$comparison\t$myComparison\n";
	if ($comparison eq $myComparison){
	    for (my $i=1;$i<=$#data;$i++){
		my $sample = $samples{$i};
		$sample =~ s/\_/\./g;
		if ($data[$i] == 1){
		    $group1samples{$sample} = 1;
		} elsif ($data[$i] == -1){
		    $group2samples{$sample} = 1;
		}
	    }
	}
    }
}

my @group1samples = keys(%group1samples);
my @group2samples = keys(%group2samples);
print STDERR "Averaging data in $myComparison for $dataTable\n";
my $group1size = $#group1samples +1;
my $group2size = $#group2samples +1;
print STDERR "Group1 has $group1size samples:\n@group1samples\n";
print STDERR "Group2 has $group2size samples:\n@group2samples\n";

my @group1indices;
my @group2indices;
open (IN,$dataTable);
while (<IN>){
    chomp $_;
    if (($_ =~ /^Read/) || ($_ =~ /^Avg/) || ($_ =~ /^Seed/)){
	@labels = split(/\t/,$_);
	for (my $i=0;$i<=$#labels;$i++){
	    my $label = $labels[$i];
	    if (exists($group1samples{$label})) {
		push (@group1indices,$i);
	    }
	    if (exists($group2samples{$label})) {
		push (@group2indices,$i);
	    }
	}
	print STDERR "Group1 columns: @group1indices\n";
	print STDERR "Group2 columns: @group2indices\n";
	print $_;
	print "\t$myComparison.Group1\t$myComparison.Group2\n";
    } else {
	my @data = split(/\t/,$_);
	my $group1sum = 0;
	my $group2sum = 0;
	foreach my $j (@group1indices){	       
	    $group1sum += $data[$j];
	}
	foreach my $k (@group2indices){
	    $group2sum += $data[$k];
	}
	$group1sum = $group1sum / $group1size;
	$group2sum = $group2sum / $group2size;
	print "$_\t$group1sum\t$group2sum\n";
#	print STDERR "$data[0]\t$group1sum\t$group2sum\n";
    }
}
