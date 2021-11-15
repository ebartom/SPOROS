#!/software/activeperl/5.16/bin/perl -w

use Getopt::Long qw(GetOptions);
use List::Util qw(max);
use List::MoreUtils qw(uniq);
use File::Basename;
use strict;
use utf8;
use warnings;
no warnings 'uninitialized';

unless (@ARGV) {
    print "\nUsage: buildSPOROSpipeline.pl\n\n-d\t\t<workingDirectory>\n-o\t\t<organism>\n-id\t\t<project id>\n-b\t\t<barcodes>\n-env\t\t<environmentalVariables>\n\n";
}

# Set up the environment for analysis.
my $questSporos="/projects/b1069/";
my $NGSbartom="/projects/p20742/tools/bin/";
my $pipeline="$questSporos/smallRNAscripts/";
my $usefulFiles="$questSporos/usefulFiles/";
my $pythonENV="/projects/p20742/envs/";
my $envVar="internal";
my $organism = "human";
my $project="";
my $workdir="";
my $comparisons="";
my $barcodes="";
my $account= "b1042";
my $queue = "genomics";
my $walltime = "48:00:00";
my $memory = 50000;
my $fastq = "";
my $table = "";
my $threads = 24;
my $minSum = 2;
my $runTrim = 1;
my $type = ""; # This is for UMId / notUMId if sample is pulldown.
my $sequenceType = "sRNA";
my $justComparisons = 0;

# Read in the command line arguments.
GetOptions('organism|o=s' => \$organism,
	   'project|id=s' => \$project,    
	   'workdir|w=s' => \$workdir,
	   'fastqdir|f=s' => \$fastq,
	   'table|t=s' => \$table,
	   'comparisons|c=s' => \$comparisons,
	   'sequencetype|st=s' => \$sequenceType,
	   'minSum|ms=s' => \$minSum,
	   'runTrim|rt=s' => \$runTrim,
	   'justComparisons|jc=s' => \$justComparisons,
	   'envVar|ev=s' => \$envVar,
	   'barcodes|b=s' => \$barcodes   ) ;


if ($workdir eq ""){
	$workdir = "$questSporos/$project\*";
}
if ($fastq eq ""){
    $fastq = "$workdir/fastq";
}

if ($envVar ne "internal"){
    print STDERR "Use $envVar shell script to set local paths for environmental variables.\n";
    open (ENV,$envVar);
    while(<ENV>){
	chomp $_;
	if ($_ =~ /^pipeline=([\w\/\-\.\_]+)$/){
	    $pipeline = $1;
	}
	if ($_ =~ /^pythonENV=([\w\/\-\.\_]+)$/){
	    $pythonENV = $1;
	}
	if ($_ =~ /^usefulFiles=([\w\/\-\.\_]+)$/){
	    $usefulFiles = $1;
	}
    }
}

if (-e $workdir){
    print STDERR "$workdir already exists.\n";
} else {
    `mkdir $workdir`;
}

my $shellscript = "$workdir/runSPOROSpipeline.$project.sh";
if ($justComparisons == 1){
    $shellscript = "$workdir/runSPOROSpipeline.$project.justComparisons.sh";
} 
    
print STDERR "Creating Shell script $shellscript\n";
open(SH,">$shellscript");

# This header is configured for Northwestern University's High
# Performance Compute Cluster, Quest and its SLURM scheduler.
# It should be adjusted for use on other clusters as necessary.
my $header = "#!/bin/bash\n";
$header .= "#SBATCH -A $account\n";
$header .= "#SBATCH -p $queue\n";
$header .= "#SBATCH -t $walltime\n";
$header .= "#SBATCH -m a\n"; # only email user if job aborts
$header .= "#SBATCH --mem=$memory\n";
$header .= "#SBATCH --chdir=$workdir\n";
$header .= "#SBATCH -o \"\%x.o\%j\"\n";
$header .= "#SBATCH --job-name=$project.SPOROS\n";
$header .= "#SBATCH --nodes=1\n";
$header .= "#SBATCH -n $threads\n";

print SH "$header\n";
print SH "pipeline=$pipeline\n";
print SH "project=$project\n";
print SH "workdir=$workdir\n";
print SH "usefulFiles=$usefulFiles\n";
print SH "pythonENV=$pythonENV\n";


if ($table eq ""){
    print SH "fastq=$fastq\n";
}
print SH "\n";
print SH "module load perl/5.16\n";
print SH "module load python/anaconda\n";
print SH "module load R/3.2.1\n";
print SH "module load bedtools/2.29.2\n";
print SH "module load blast/2.7.1\n";
print SH "\n";
print SH "cd \$workdir\n";

print SH "# Starting analysis for SPOROS pipeline.\n";

if ($justComparisons == 1){
    print SH "# Assuming that the preliminary analysis and totalCounts is done and those directories correctly populated.  Generating only the differential analysis.\n";
} elsif ($justComparisons == 0){
    print SH "date\n\n";
    print SH "mkdir \$workdir/prelimAnalysis\n";
    if ($table eq ""){
	print SH "mkdir \$workdir/\$project.fastq\n";
	print SH "\n";
    } elsif ($table ne ""){
	print SH "cd \$workdir/prelimAnalysis\n";
	print SH "cp $table \$workdir/prelimAnalysis\n";
    }
    
    if ($barcodes ne ""){
	my @barcodes = split(/\,/,$barcodes);
	print SH "\n# Barcode files for demultiplexing found.  @barcodes\n";
	print SH "\n# Each barcode needs to match a fastq file in the $fastq directory, with the same file prefix.\n";
	print SH "mkdir \$workdir/demultiplex\n";
	print SH "cd \$workdir/demultiplex\n";
	my $sampleGroup = "";
	for my $barcode (@barcodes){
	    print SH "\ncp $barcode \$workdir/demultiplex/ \n";
	    if ($barcode =~ /\/?([\w\_\-\.]+).barcodes.txt/){
		$sampleGroup = $1;
	    print SH "# Demultiplexing samples for sample group $sampleGroup\n";
		print SH "perl $NGSbartom/trim_galore \$fastq/$sampleGroup.fastq.gz --length 12 --dont_gzip -o \$workdir/demultiplex/ \n";
		print SH "perl \$pipeline/splitFastQwithTable.pl $sampleGroup\_trimmed.fq $sampleGroup.barcodes.txt 0\n";
		print SH "mv *.fastq \$workdir/\$project.fastq/\n";
		print SH "gzip \$workdir/\$project.fastq/*.fastq \&\n";
		print SH "gzip \$workdir/demultiplex/$sampleGroup*.fq \&\n";
	    }
	}
	print SH "wait\n";
	print SH "\n# Finished de-multiplexing.\n";
	print SH "date\n";
	print SH "cd \$workdir\n";
    } elsif (($barcodes eq "") && ($table eq "")) {
	if ($runTrim == 1){
	    print SH "\n\# Trim the reads to remove any standard illumina adaptors and remove sequences < 18 bp or > 25 bp.\n";
	    print SH "cd $fastq\n";
	    print SH "for f in \*.fastq.gz\n";
	    print SH "do\n";
	    print SH "\techo \$f\n";
	    print SH "\tperl $NGSbartom/trim_galore \$f --length 18 -o \$workdir/\$project.fastq/ --dont_gzip --max_length 25 \&\n";
	    print SH "done\n\n";
	    print SH "\# Wait for trimming to finish.\n";
	    print SH "wait\n";
	    print SH "cp \$fastq/*.fastq.gz \$workdir/\$project.fastq/\n";
	    print SH "\n# Finished trimming.\n";
	    print SH "date\n";
	    print SH "cd \$workdir\n";
	} elsif ($runTrim == 0){
	    print SH "\n# Copy reads over from input fastq folder to project fastq.\n";
	    print SH "cp \$fastq/*.fastq.gz \$workdir/\$project.fastq/\n";
	}
    }
    if ($table eq ""){
	print SH "\n";
	print SH "\# Building counts table from raw reads.\n";
	print SH "mkdir \$workdir/prelimAnalysis\n";
	print SH "cd \$workdir/prelimAnalysis\n";
	print SH "\# Pull out the individual reads for read-based analyses.\n";
	print SH "for f in \$workdir/$project.fastq/\*.f*q.gz\n";
	print SH "do\n";
	print SH "\t\# Determine the sample name based on the fastq file name.\n";
	print SH "\techo \$f\n";
	print SH "\tfile=\$\(basename \$f\)\n";
	print SH "\techo \$file\n";
	print SH "\tsample=\$\{file\%\%.gz\}\n";
	print SH "\tsample=\$\{sample\%\%\_trimmed.fq\}\n";
	print SH "\tsample=\$\{sample\%\%.fastq\}\n";
	print SH "\n\t# Print the sample name for the fastq file.\n";
	print SH "\techo \$sample\n";
	print SH "\tgunzip \$f\n";
	print SH "\n\t# Pull out just the read sequences, leaving quality scores behind.\n";
	if ($barcodes eq ""){
	    if ($runTrim == 1){
		print SH "\tgrep \^\@ -A 1 \$workdir/\$project.fastq/\$sample\\_trimmed.fq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > \$workdir/prelimAnalysis/\$sample.justReads.fastq\n";
	    }elsif ($runTrim == 0){
		print SH "\tgrep \^\@ -A 1 \$workdir/\$project.fastq/\$sample.fastq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > \$workdir/prelimAnalysis/\$sample.justReads.fastq\n";
	    }
	} elsif ($barcodes ne "") {
	    print SH "\tgrep \^\@ -A 1 \$workdir/\$project.fastq/\$sample.fastq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > \$workdir/prelimAnalysis/\$sample.justReads.fastq\n";
	}
	print SH "\n\t# Combine duplicate reads and add read counts to each line.\n";
	print SH "\tsort \$workdir/prelimAnalysis\/\$sample.justReads.fastq | uniq -c | sort -nr > \$workdir/prelimAnalysis\/\$sample.justReads.uniqCounts.txt\n";
	print SH "\tgzip \$workdir/\$project.fastq/\$sample*q &\n";
	print SH "\tgzip \$workdir/prelimAnalysis/\$sample.justReads.fastq &\n";
	print SH "done\n";
	print SH "\n";
	print SH "cd \$workdir/prelimAnalysis\n";
	if ($barcodes ne ""){
	    print SH "\n# We assume that samples with barcodes also have UMI sequences of the pattern 5'-NNNNreadNN-3'.  If this is not the case, fastq files should be demultiplexed before starting SPOROS, and it should be run without demultiplexing barcodes.\n"; 
	    $type = "notUMId.";
	    print SH "\n# Process UMI sequences to split out reads into two sets, UMId and notUMId.\n";
	    print SH "for f in \*.justReads.uniqCounts.txt\n";
	    print SH "do\n";
	    print SH "\t\# Extract the sample name from the file name.\n";
	    print SH "\tfile=\$\(basename \$f\)\n";
	    print SH "\techo \$file\n";
	    print SH "\tsample=\$\{file\%\%.justReads.uniqCounts.txt\}\n";
	    print SH "\techo \$sample\n";
	    print SH "\n\t\# Generate two files, and remove UMI sequences from reads.\n";
	    print SH "\t# notUMId: A six nucleotide Unique Molecular Identifier (UMI) is removed from each read, and the numbers of reads associated with each UMI are summed to generate the raw read count for a particular sample.\n";
	    print SH "\t# For the UMId files, the number of unique UMIs is counted, regardless of how frequently each UMI is seen.\n";
	    print SH "\tperl \$pipeline\/processUMIs.pl \$sample.justReads.uniqCounts.txt\n"; 
	    print SH "done\n";
	}
	print SH "\n# Remove adapter sequences from raw read counts.\n";
    	print SH "for f in \*.justReads.uniqCounts*txt\n";	
	print SH "do\n";
	print SH "\t# Adapter reads are identified as containing the substring TCCGACGATC followed by 3-5 random nucleotides and removed from analysis.\n";
	print SH "\tperl \$pipeline\/deleteAdapterReads.pl \$f > \$f.noAdapter.txt\n";
	print SH "done\n";
	print SH "\n# Make tables from the reads.\n";
	print SH "# Counts are tabulated for each read and sample. The rawCounts table has raw un-normalized data.\n";
	print SH "perl \$pipeline/combineCountsNotNorm.pl justReads.uniqCounts.".$type."txt.noAdapter.txt 1 > \$project.noAdapter.".$type."rawCounts.table.txt \n";
#	print SH "# At this step, remove reads / rows with a count sum < the number of samples in the table. The last column in the table will be the row total.\n";
#	print SH "awk  \' \$NF \>= \(NF-2\) \' \$project.noAdapter.".$type."rawCounts.table.txt > \$project.noAdapter.".$type."rawCounts.minSumN.table.txt\n";
	$table = "$project.noAdapter.".$type."rawCounts.table.txt";
    }

    print SH "\n# Using table of raw counts $table for the rest of the pipeline.\n";
}
# If not running the full pipeline, and critical variables are empty, define them.
if ($table eq ""){
    if ($type eq ""){
	if ($barcodes ne ""){
	    $type = "notUMId.";
	}
    }
    $table = "$project.noAdapter.".$type."rawCounts.table.txt";
}
my $tablePrefix = "";
if (basename($table) =~ /^([\w\_\.\-]+).txt$/){
    $tablePrefix = $1;
}
print SH "cd \$workdir/prelimAnalysis\n";
print SH "tablePrefix=$tablePrefix\n";

if ($justComparisons == 0){
    print SH "# Remove reads / rows with a count sum < the number of samples in the table. The last column in the table should be the row total.\n";
    print SH "awk  \' \$NF \>= \(NF-2\) \' $table > \$tablePrefix.minSumN.txt\n";
    print SH "\n# Make a fasta formatted file from the reads in the table\n";
    print SH "awk \' \$1 \!\~ \"Read\" \{printf \"\>\%s.\%d\\n\%s\\n\",\$1,\$NF,\$1\}\' $tablePrefix.minSumN.txt > \$tablePrefix.reads.fa\n";
    print SH "\n# Reads are blasted (blastn-short) against custom blast databases created from all processed miRNAs and from the most recent RNA world databases, using the sequences appropriate for the organism (human or mouse). \n";
    if ($organism eq "human"){
	print SH "# Blast reads against human miRNAs\n";
	print SH "blastn -task blastn-short -query \$tablePrefix.reads.fa -db $usefulFiles/humanMiRNAs.nr -outfmt 7 -num_threads $threads > \$tablePrefix.miRNAs.nr.$organism.blast.txt \& \n";
	print SH "\# Blast reads against human RNA world.\n";
	print SH "blastn -task blastn-short -query \$tablePrefix.reads.fa -db $usefulFiles/human_and_virus_vbrc_all_cali_annoDB_sep2015.proc -outfmt 7 -num_threads $threads > \$tablePrefix.RNAworld.$organism.blast.txt &\n";
    } elsif ($organism eq "mouse"){
	print SH "\# Blast reads against mouse miRNAs\n";
	print SH "blastn -task blastn-short -query \$tablePrefix.reads.fa -db $usefulFiles/mouseMiRNAs.nr -outfmt 7 -num_threads $threads > \$tablePrefix.miRNAs.nr.$organism.blast.txt \& \n";
	print SH "\# Blast reads against mouse RNA world.\n";
	print SH "blastn -task blastn-short -query \$tablePrefix.reads.fa -db $usefulFiles/mouse_annoDB_piRNA_all_oct2013.proc -outfmt 7 -num_threads $threads > \$tablePrefix.RNAworld.$organism.blast.txt &\n";
    }
print SH "wait\n";
    print SH "\n\# After blast runs are complete, resume analysis.\n";
    print SH "\n\# Analyze blast data.\n";
    print SH "# Filter blast results\n";
    print SH "# Processed miRNA hits were filtered for 100% identity, and a match length of at least 18 bp.\n";
    print SH "awk \' \$3 == 100 && \$8 > \$7 && \$10 > \$9 && \$4 >= 18 && \$9 < 10 \' \$tablePrefix.miRNAs.nr.$organism.blast.txt  > \$tablePrefix.miRNAs.nr.blast.filtered.txt\n";
    print SH "# RNAworld hits were filtered for at least 95% identity, where the match is within the first 9 bp of the sequencing reads.\n";
    print SH "awk \' \$3 >= 95 && \$8 > \$7 && \$10 > \$9 && \$7 < 10\'  \$tablePrefix.RNAworld.$organism.blast.txt  > \$tablePrefix.RNAworld.blast.filtered.txt\n";
    print SH "\n# Remove reads from table that have hits against RNAworld adapter sequences.\n";
    print SH "perl \$pipeline/removeRNAworldAdapterSeq.pl \$tablePrefix.minSumN.txt \$tablePrefix.RNAworld.blast.filtered.txt > \$tablePrefix.RNAworldClean.txt\n";
    print SH "cp \$tablePrefix.RNAworldClean.txt $table\n";
    print SH "\n# Normalize raw read counts to counts per million, removing reads with fewer than $minSum counts across all samples. (Note that reads with fewer reads than the number of samples have already been removed.)\n";
    print SH "perl \$pipeline/normalizeRawCountsToCPM.pl \$tablePrefix.RNAworldClean.txt $minSum > \$tablePrefix.normCounts.txt\n";
    
    print SH "\n # Add seeds, toxicities, blast hits to normalized counts table.\n";
    print SH "perl $usefulFiles/addAllToxicities.pl $usefulFiles/speciesToxes.txt \$tablePrefix.normCounts.txt > \$tablePrefix.normCounts.withTox.txt\n";
    print SH "perl \$pipeline/addBLASTresults.pl \$tablePrefix.normCounts.withTox.txt \$tablePrefix.miRNAs.nr.blast.filtered.txt > \$tablePrefix.normCounts.withTox.withMiRNAblast.txt\n";
    print SH "perl \$pipeline/addBLASTresults.pl \$tablePrefix.normCounts.withTox.withMiRNAblast.txt \$tablePrefix.RNAworld.blast.filtered.txt > \$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt\n";
    print SH "perl \$pipeline/truncateBLASTresults.pl \$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt > \$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt\n";
    print SH "\n # Add seeds, toxicities, blast hits to raw counts table.\n";
    print SH "perl $usefulFiles/addAllToxicities.pl $usefulFiles/speciesToxes.txt \$tablePrefix.RNAworldClean.txt > \$tablePrefix.withTox.txt\n";
    print SH "perl \$pipeline/addBLASTresults.pl \$tablePrefix.withTox.txt \$tablePrefix.miRNAs.nr.blast.filtered.txt > \$tablePrefix.withTox.withMiRNAblast.txt\n";
    print SH "perl \$pipeline/addBLASTresults.pl \$tablePrefix.withTox.withMiRNAblast.txt \$tablePrefix.RNAworld.blast.filtered.txt > \$tablePrefix.withTox.withMiRNAandRNAworld.blast.txt\n";
    print SH "perl \$pipeline/truncateBLASTresults.pl \$tablePrefix.withTox.withMiRNAandRNAworld.blast.txt > \$tablePrefix.withTox.withMiRNAandRNAworld.blast.trunc.txt\n";
    
#    print SH "\n# Pull out smaller subsets of these tables for easy manipulation in Excel.\n";
    #print SH "perl \$pipeline/thresholdNormTotal.pl \$tablePrefix.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.txt 20 > \$tablePrefix.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.minSum20.txt\n";
#    print SH "perl \$pipeline/thresholdNormTotal.pl \$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt 6 > \$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.minSum6.txt\n";
    print SH "\n# This \"normCounts\" file is used as the basis for many further analyses in the Peter Lab, and is the first of the official output files for the SPOROS pipeline.  For simplicity, we will rename it A_normCounts.\$project.txt and put it in a new directory with the other output files.\n";
    print SH "mkdir \$workdir/totalCounts\n";
    #print SH "cp \$workdir/prelimAnalysis/\$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt \$workdir/totalCounts/A_normCounts.\$project.txt\n";
    print SH "perl -pe \"s/$tablePrefix.RNAworld.blast.filtered.txt/RNAworld/g\"  \$workdir/prelimAnalysis/\$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe \"s/$tablePrefix.miRNAs.nr.blast.filtered.txt/miRNA/g\" > \$workdir/totalCounts/A_normCounts.\$project.txt\n";
    print SH "# Also copy over rawCounts version.\n";
    print SH "perl -pe \"s/$tablePrefix.RNAworld.blast.filtered.txt/RNAworld/g\"  \$workdir/prelimAnalysis/\$tablePrefix.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe \"s/$tablePrefix.miRNAs.nr.blast.filtered.txt/miRNA/g\" > \$workdir/totalCounts/A_rawCounts.\$project.txt\n";
    print SH "cd \$workdir/totalCounts\n";
    print SH "\n# Collapse counts for seed sequences from the same RNA species.\n";
    print SH "perl \$pipeline/collapseSpecies.pl A_normCounts.\$project.txt > B_collapsed.\$project.txt\n";
    print SH "\n# Split the analysis by sequence type and generate a directory for each and run separately.\n";
    print SH "for seqtype in sRNA miRNA\n";
    print SH "do\n";
    print SH "\tmkdir \$workdir/totalCounts/\$seqtype\n";
    print SH "\tcd \$workdir/totalCounts/\$seqtype\n";
    print SH "\n\t# Collapse counts for seed sequences, regardless of RNA species of origin, for all \$seqtype.\n";
    print SH "\tperl \$pipeline/collapseToxicityBins.pl ../B_collapsed.\$project.txt $organism \$seqtype\n";
    print SH "\t# This generates files C_binned.\$seqtype.$organism.\$project.txt and intermediary file Int_seedKeyed.\$seqtype.$organism.\$project.txt.\n";
    print SH "\n\t# Expand tox counts for each sample (and the average of any replicates) into lines in a text file, for boxplots.\n";
    print SH "\tperl \$pipeline/expandToxesFromSeedCollapsed.pl ../B_collapsed.\$project.txt $organism \$seqtype \$project\n";
    print SH "\t# This generates files that start with the word D_toxAnalysis, one for each sample and average of replicate samples.\n";
    print SH "\n\t# Expand seed counts for each samples (and the average of any replicates) into lines in a text file, for weblogo.\n";
    print SH "\tperl \$pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.\$seqtype.$organism.\$project.txt \$seqtype \$project\n";
    print SH "\t# This generates files that start with the words E_seedAnalysis and F_seedExpand, one for each sample and average of replicate samples.\n";
    print SH "\n\t# Run weblogo.\n";
    print SH "\tmodule load python/anaconda3.6\n";
    print SH "\tsource activate $pythonENV/weblogo-py38\n";
    print SH "\tfor f in E_seedAnalysis*.txt\n\tdo\n";
    print SH "\t\techo \$f\n";
    print SH "\t\tsample=\$\{f\%\%.txt\}\n";
#    print SH "\t\tweblogo -f \$f -A rna --show-xaxis NO --show-yaxis NO -U probability -F pdf -o \$sample.pdf \\\n\t\t\t--color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t\t\t--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1\n";
    print SH "\t\tweblogo -f \$f -A rna --show-xaxis NO --show-yaxis NO -U probability -F png_print -o \$sample.png \\\n\t\t\t--color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t\t\t--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1\n";
    print SH "\t\tdone \n";
    print SH "\tsource deactivate\n\n";
    print SH "\t# Compile the toxAnalysis files into one multi-column summary.\n";
    print SH "\t# If there are averaged files, use those.\n";
    print SH "\tif ls | grep -q \"avg\"; then\n";
    print SH "\t\tfor d in D_toxAnalysis.*.avg.*.txt\n";
    print SH "\t\tdo\n";
    print SH "\t\t\techo \$d\n";
    print SH "\t\t\tsample=\${d\#\#D_toxAnalysis.\}\n";
    print SH "\t\t\tsample=\${sample\%\%.human.$project.txt}\n";
    print SH "\t\t\techo \$sample | cat > test.\$sample.tox\n";
    print SH "\t\t\tcat \$d >> test.\$sample.tox\n";
    print SH "\t\tdone\n";
    print SH "\t# If there are no averaged files, just use all of the individual ones.\n";
    print SH "\telse\n";
    print SH "\t\tfor d in D_toxAnalysis.*.txt\n";
    print SH "\t\tdo\n";
    print SH "\t\t\techo \$d\n";
    print SH "\t\t\tsample=\${d\#\#D_toxAnalysis.\}\n";
    print SH "\t\t\tsample=\${sample\%\%.human.$project.txt}\n";
    print SH "\t\t\techo \$sample | cat > test.\$sample.tox\n";
    print SH "\t\t\tcat \$d >> test.\$sample.tox\n";
    print SH "\t\tdone\n";
    print SH "\tfi\n";
    print SH "\tpaste test.\*.tox > D_toxAnalysis.combined.\$seqtype.txt\n";
    print SH "\trm test.\*.tox\n\n";
    print SH "\t# Plot a very basic box plot summarizing the tox distributions in the different samples.\n";
    print SH "\tRscript $pipeline/plotToxBoxPlot.R --toxFile=D_toxAnalysis.combined.\$seqtype.txt\n\n";
    print SH "\t# Analysis done for \$seqtype.\n";
    print SH "\tdone\n";
    print SH "date\n";
    print SH "\n# Analysis done for total counts.\n";
}

if ($comparisons ne ""){
    print SH "\n# Running differential read analysis.\n";
    print SH "\nmkdir \$workdir/differential\n";
    print SH "module load R/3.2.1\n";
    print SH "cd \$workdir/prelimAnalysis\n";
    print SH "\n# Run EdgeR on the comparisons file to identify differentially abundant reads.\n";
    print SH "# (Note that assembly is not actually used here, because these are reads, not gene IDs)\n";
    print SH "perl -pe \"s\/.justReads.uniqCounts.".$type."txt\/\/g\" $table | perl -pe \"s/.noAdapter.txt//g\" | perl -pe \"s\/\_S\\d\+\_R1\_001\/\/g\" | perl -pe \"s/.geo.txt//g\" > \$tablePrefix.relabeled.txt\n";
    print SH "Rscript \$pipeline/runEdgeRrnaSeq.reads.R --countFile=\$tablePrefix.relabeled.txt --numCores=$threads --runMDS=1 --filterOff=0 --assembly=hg38.mp\n";
    print SH "Rscript \$pipeline/runEdgeRrnaSeq.reads.R --countFile=\$tablePrefix.relabeled.txt --numCores=$threads --runMDS=0 --filterOff=0 --assembly=hg38.mp --comparisonFile=$comparisons\n";
    print SH "\n# Annotate each EdgeR result file with toxicities and blast results.\n";
    print SH "for edgeR in *relabeled.txt.edgeR.txt\n";
    print SH "do\n";
    print SH "\techo \$edgeR\n";
    print SH "\tcomp=\$\{edgeR\%\%.$tablePrefix.relabeled.txt.edgeR.txt\}\n";
    print SH "\techo \$comp\n";
    print SH "\tcd \$workdir/prelimAnalysis/\n";
    print SH "\tperl -pe \"s\/\\\"\/\/g\" \$edgeR > \$comp.txt\n";
    print SH "\n\t# Add seed toxicity and source sequences to edgeR results.\n";
    print SH "\tperl $usefulFiles/addAllToxicities.pl $usefulFiles/speciesToxes.txt \$comp.txt > \$comp.withTox.txt\n";
    print SH "\tperl \$pipeline/addBLASTresults.pl \$comp.withTox.txt \$tablePrefix.miRNAs.nr.blast.filtered.txt > \$comp.withTox.withMiRNAblast.txt\n";
    print SH "\tperl \$pipeline/addBLASTresults.pl \$comp.withTox.withMiRNAblast.txt \$tablePrefix.RNAworld.blast.filtered.txt > \$comp.withTox.withMiRNAandRNAworld.blast.txt\n";
    print SH "\tperl \$pipeline/truncateBLASTresults.pl \$comp.withTox.withMiRNAandRNAworld.blast.txt > \$comp.withTox.withMiRNAandRNAworld.blast.trunc.txt\n";
    print SH "\n\t# Rename file for simplicity and get organized.\n";
    print SH "\tmkdir \$workdir/differential/\$comp\n";
    #print SH "\tcp \$workdir/prelimAnalysis/\$comp.withTox.withMiRNAandRNAworld.blast.trunc.txt \$workdir/differential/\$comp/A_normCounts.\$comp.all.txt\n";
    print SH "\tperl -pe \"s/$tablePrefix.RNAworld.blast.filtered.txt/RNAworld/g\"  \$workdir/prelimAnalysis/\$comp.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe \"s/$tablePrefix.miRNAs.nr.blast.filtered.txt/miRNA/g\" > \$workdir/differential/\$comp/A_normCounts.\$comp.all.txt\n";
    print SH "\tcd \$workdir/differential/\$comp\n";
    print SH "\t# Create a subdirectory for both an adjusted p-value comparison and a p-value comparison, and then run both of them within the relevant sub-directory\n";
    print SH "\tfor stat in adjp pvalue\n";
    print SH "\tdo\n";
    print SH "\t\tmkdir \$workdir/differential/\$comp/\$stat\n";	
    print SH "\n\t\t# Separate out reads that are differentially regulated at a \$stat <= 0.05 and a absolute logFC > 0.585\n";
    print SH "\t\tcd \$workdir/differential/\$comp/\$stat\n";
    print SH "\t\tperl \$pipeline/thresholdDifferential2.pl ../A_normCounts.\$comp.all.txt \$stat 0.585 0.05\n";
    print SH "\t\t# This generates files A_normCounts.\$comp.diff.txt with all differentially regulated reads.\n\n";
    print SH "\n\t\t# Collapse counts for seed sequences from the same RNA species.\n";
    print SH "\t\tperl \$pipeline/collapseSpecies.pl A_normCounts.\$comp.diff.txt > B_collapsed.\$comp.diff.txt\n";
    print SH "\n\t\t# Guided by comparisons file, find the average counts across the two groups indicated for the comparison \$comp, for differential reads\n";
    print SH "\t\tperl \$pipeline/averageWithinGroups.pl $comparisons B_collapsed.\$comp.diff.txt \$comp > B_collapsed.\$comp.diff.avgd.txt\n";
    print SH "\t\t# Split downstream analysis depending on whether the reads from miRNAs or sRNAs.\n";
    print SH "\t\tfor seqtype in sRNA miRNA\n";
    print SH "\t\tdo\n";
    print SH "\t\t\tmkdir \$workdir/differential/\$comp/\$stat/\$seqtype\n";
    print SH "\t\t\tcd \$workdir/differential/\$comp/\$stat/\$seqtype\n";
    print SH "\n\t\t\t# Collapse counts for differential and up-regulated read sequences, regardless of RNA species of origin, as long as they are \$seqtype.\n";
    print SH "\t\t\tperl \$pipeline/collapseToxicityBins.pl ../B_collapsed.\$comp.diff.avgd.txt $organism \$seqtype\n";
    print SH "\t\t\t# This generates files C_binned.\$seqtype.$organism.\$comp.diff.avgd.txt and intermediary file Int_seedKeyed.\$seqtype.$organism.\$comp.diff.avgd.txt.\n\n";
    print SH "\t\t\t# Subtract read counts in Group1 from Group2 to find the difference between the two.\n";
    print SH "\t\t\tperl \$pipeline/deltaBetweenGroups.pl Int_seedKeyed.\$seqtype.$organism.\$comp.diff.avgd.txt\n";
    print SH "\t\t\t# This generates files Int_seedKeyed.\$seqtype.$organism.\$comp.delta.up.txt and Int_seedKeyed.\$seqtype.$organism.\$comp.delta.up.txt \n";
    print SH "\n\t\t\t# Expand delta seed counts for seeds that went either up or down into lines in a text file, for weblogo.\n";
    print SH "\t\t\tperl \$pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.\$seqtype.$organism.\$comp.delta.up.txt \$seqtype up\n";
    print SH "\t\t\tperl \$pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.\$seqtype.$organism.\$comp.delta.dn.txt \$seqtype dn\n";
    print SH "\n\t\t\t# Expand tox counts for each delta up or delta down file into lines in a text file, for boxplots.\n";
    print SH "\t\t\tperl \$pipeline/expandToxesFromSeedKeyed.pl Int_seedKeyed.\$seqtype.$organism.\$comp.delta.up.txt up\n";
    print SH "\t\t\tperl \$pipeline/expandToxesFromSeedKeyed.pl Int_seedKeyed.\$seqtype.$organism.\$comp.delta.dn.txt dn\n";
    #print SH "\t\t\tperl \$pipeline/expandToxesFromSeedKeyed.pl B_collapsed.\$comp.delta.dn.txt dn\n";
##    print SH "\t\tperl \$pipeline/expandToxesFromSeedCollapsed.pl B_collapsed.\$comp.diff.avgd.txt $organism \$seqtype diff\n";

    print SH "\n\t\t\t# Run weblogo.\n";
    print SH "\t\t\tmodule load python/anaconda3.6\n";
    print SH "\t\t\tsource activate $pythonENV/weblogo-py38\n";
    print SH "\t\t\tfor f in E_seedAnalysis*Delta*txt\n";
    print SH "\t\t\tdo\n";
    print SH "\t\t\t\techo \$f\n";
    print SH "\t\t\t\tsample=\$\{f\%\%.txt\}\n";
#    print SH "\t\t\t\tweblogo -f \$sample.txt -A rna --show-xaxis NO --show-yaxis NO -U probability -F pdf -o \$sample.pdf \\\n\t\t\t\t\t--color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t\t\t\t\t--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1\n";
print SH "\t\t\t\tweblogo -f \$sample.txt -A rna --show-xaxis NO --show-yaxis NO -U probability -F png_print -o \$sample.png \\\n\t\t\t\t\t--color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t\t\t\t\t--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1\n";
    print SH "\t\t\t\tdone\n";
    print SH "\t\t\tsource deactivate\n\n";
    print SH "\t\t\t# Compile the delta toxAnalysis files into one multi-column summary.\n";
    print SH "\t\t\tfor d in D_toxAnalysis.delta.*.txt\n";
    print SH "\t\t\tdo\n";
    print SH "\t\t\t\techo \$d\n";
    print SH "\t\t\t\tsample=\${d\#\#D_toxAnalysis.\}\n";
    print SH "\t\t\t\tsample=\${sample\%\%.human.$project.txt}\n";
    print SH "\t\t\t\techo \$sample | cat > test.\$sample.tox\n";
    print SH "\t\t\t\tcat \$d >> test.\$sample.tox\n";
    print SH "\t\t\tdone\n";
    print SH "\t\t\tpaste test.\*.tox > D_toxAnalysis.combined.\$seqtype.txt\n";
    print SH "\t\t\trm test.\*.tox\n\n";
    print SH "\t\t\t# Plot a very basic box plot summarizing the tox distributions in the different samples.\n";
    print SH "\t\t\tRscript $pipeline/plotToxBoxPlot.R --toxFile=D_toxAnalysis.combined.\$seqtype.txt\n";
    print SH "\t\t\tdone\n";
    print SH "\t\tdone\n";
    print SH "\tdone\n";
    print SH "done\n";
}
