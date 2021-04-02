#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use List::Util qw( min max );

my $fastQfile = $ARGV[0];
my $barcodes = $ARGV[1];
my $spacer = $ARGV[2];

my $minLen = 6;

if ($#ARGV <0){
    die "Usage:  $0 <fastQfile> <barcodeTable> [spacer]\n";
}
if (!($spacer)){$spacer = 0;}
print STDERR "Spacer: $spacer\n";

my %samples;
open(IN,$barcodes);
while(<IN>){
    chomp $_;
    my ($barcode,$sample) = split(/\t/,$_);
    if ($barcode ne "barcode"){
	print STDERR "\"$barcode\"\t$sample\n";
	$samples{$barcode} = $sample;
    }
}

print STDERR "Sequence file: $fastQfile\n";
my $fastqlabel = $fastQfile;
$fastqlabel =~ s/.fastq$//g;
$fastqlabel =~ s/.fq$//g;
foreach my $barcode (keys(%samples)){
    my $sample = $samples{$barcode};
    open(FQ,$fastQfile);
    open (OUT1,">$sample.fastq");
    my $lineCount = 5;
    my ($header,$seq,$qual,$header2);
    my ($catchStop,$catchQual,$catchStart);
    my $readCount = 0;
    while(<FQ>){
	chomp $_;
	#  print STDERR "$lineCount\t$_\n";
	if (($_ =~ /^\@/)&& ($lineCount == 5)){ $header =$_; $lineCount = 1;}
	elsif (($lineCount == 2)){ $seq = $_;}
	elsif (($lineCount == 3)){ $header2 = $_}
	elsif ($lineCount == 4){ 
	    $qual = $_;
	    #	    print STDERR "$qual\n";
	    #    print STDERR "$header $seq\n";	       
	    my $catch = "";
	    if ($seq =~ /([\w\.]+)$barcode$/){
		$catchStop = $+[1];
		if ($spacer > 0){
		    $catchStop -= $spacer;
		}
		$catchQual = substr($qual,0,$catchStop);
		$catch = substr($seq,0,$catchStop);
		my $catchLen = length($catch);
		my $numN = () = $catch =~ /N/gi;
		if (($catchLen >= $minLen) && ($numN < ($catchLen/2))){
		    print OUT1 "$header\n$catch\n+\n$catchQual\n";
		    $readCount++;
		} 
	    } 
	}
	#print STDERR "$lineCount\n";
	$lineCount++;
    }
    close(FQ);
    print STDERR "$sample\t$barcode\t$readCount\n";
}
close(OUT1);


sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
    
sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
