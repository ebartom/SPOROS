#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < comparisonsFile > [project] \n";
  exit;
}

my $compfile = $ARGV[0];
my $project = $ARGV[1];
if ($project){ $project = ".$project";} else {$project = "";}

my @labels;
my @data;
my %samples;

open (IN,$compfile);
while (<IN>){
    chomp $_;
#    print "line: $_\n";
    if ($_ =~ /^Comparison/){
	@labels = split(/\,/,$_);
	for (my $i=0;$i<=$#labels;$i++){
	    $samples{$i} = $labels[$i];
#	    print "$labels[$i]\n";
	}
    } else {
	@data = split(/\,/,$_);
	my $comparison = $data[0];
	my @group1samples;
	my @group2samples;
	for (my $i=1;$i<=$#data;$i++){
	    my $sample = $samples{$i};
	    $sample =~ s/\_/\./g;
	    if ($data[$i] == 1){
		push (@group1samples,$sample);
	    } elsif ($data[$i] == -1){
		push (@group2samples,$sample);
	    }
	}
	print "$comparison\n@group1samples\n@group2samples\n";
	foreach my $direction ("up","dn") {
	    my $outfile1 = "seedAnalysis.$comparison".$project.".group1.$direction.txt";
	    my $outfile2 = "seedAnalysis.$comparison".$project.".group2.$direction.txt";
	    system("rm $outfile1");
	    system("touch $outfile1");
	    system("rm $outfile2");
	    system("touch $outfile2");
	    foreach my $sample (@group1samples){
		my $syscommand = "cat seedAnalysis.$sample*$comparison*$direction.1000.txt >> $outfile1";
		print "$syscommand\n";
		system($syscommand);
	    }	    
	    foreach my $sample (@group2samples){
		my $syscommand = "cat seedAnalysis.$sample*$comparison*$direction.1000.txt >> $outfile2";
		print "$syscommand\n";
		system($syscommand);
	    }
	}
    }
}
