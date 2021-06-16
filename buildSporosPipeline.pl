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
    print "\nUsage: buildSporosScripts.pl\n\n-d\t\t<workingDirectory>\n-o\t\t<organism>\n-id\t\t<project id>\n-b\t\t<barcodes>\n\n";
}

# Set up the environment for analysis.
my $questSporos="/projects/b1069/";
my $NGSbartom="/projects/p20742/tools/bin/";
my $pipeline="$questSporos/smallRNAscripts/";
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
	   'barcodes|b=s' => \$barcodes   ) ;

if ($workdir eq ""){
	$workdir = "$questSporos/$project\*";
}
if ($fastq eq ""){
    $fastq = "$workdir/fastq";
}

if (-e $workdir){
    print STDERR "$workdir already exists.\n";
} else {
    `mkdir $workdir`;
}

my $shellscript = "$workdir/runSporosPipeline.$project.sh";
print STDERR "Creating Shell script $shellscript\n";
open(SH,">$shellscript");

my $header = "#!/bin/bash\n";
$header .= "#SBATCH -A $account\n";
$header .= "#SBATCH -p $queue\n";
$header .= "#SBATCH -t $walltime\n";
$header .= "#SBATCH -m a\n"; # only email user if job aborts
$header .= "#SBATCH --mem=$memory\n";
$header .= "#SBATCH --chdir=$workdir\n";
$header .= "#SBATCH -o \"\%x.o\%j\"\n";
$header .= "#SBATCH --job-name=$project.Sporos\n";
$header .= "#SBATCH --nodes=1\n";
$header .= "#SBATCH -n $threads\n";

print SH "$header\n";
print SH "pipeline=$pipeline\n";
print SH "project=$project\n";
print SH "workdir=$workdir\n";
if ($table eq ""){
    print SH "fastq=$fastq\n";
}
print SH "currentDate=\$\(date\ \+\%F\)\n";
print SH "\n";
print SH "module load perl/5.16\n";
print SH "module load python/anaconda\n";
print SH "module load R/3.2.1\n";
print SH "module load bedtools/2.29.2\n";
print SH "module load blast/2.7.1\n";
print SH "\n";
print SH "cd \$workdir\n";
print SH "mkdir \$workdir/readBased\n";

if ($table eq ""){
    print SH "mkdir \$workdir/\$project.fastq\n";
    print SH "\n";
} elsif ($table ne ""){
    print SH "cd \$workdir/readBased\n";
    print SH "cp $table \$workdir/readBased\n";
}

if ($barcodes ne ""){
    my @barcodes = split(/\,/,$barcodes);
    print SH "\n# Barcode files for demultiplexing found.  @barcodes\n";
    my $sampleGroup = "";
    for my $barcode (@barcodes){
	if ($barcode =~ /\/?([\w\_\-\.]+).barcodes.txt/){
	    $sampleGroup = $1;
	    print SH "\n\# Demultiplexing samples for sample group $sampleGroup\n";
	    print SH "perl $NGSbartom/trim_galore \$fastq/$sampleGroup.fastq.gz --length 12 --dont_gzip\n";
	    print SH "perl \$pipeline/splitFastQwithTable.pl $sampleGroup\_trimmed.fq $sampleGroup.barcodes.txt 0\n";
	    print SH "mv *.fastq \$workdir/\$project.fastq/\n";
	    print SH "gzip \$workdir/\$project.fastq/*.fastq \&\n";
	}
    }
    print SH "wait\n";
} elsif (($barcodes eq "") && ($table eq "")) {
    if ($runTrim == 1){
	print SH "\n\# Trim the reads to remove any standard illumina adaptors.\n";
	print SH "cd $fastq\n";
	print SH "for f in \*.fastq.gz\n";
	print SH "do\n";
	print SH "\techo \$f\n";
	print SH "\tperl $NGSbartom/trim_galore \$f --length 16 -o \$workdir/\$project.fastq/ --dont_gzip \&\n";
	print SH "done\n\n";
	print SH "\# Wait for trimming to finish.\n";
	print SH "wait\n";
#	print SH "\# Move trimmed files into \$project.fastq directory.\n";
#	print SH "mv \*trimmed.fq.gz \$workdir/\$project.fastq/\n";
    } elsif ($runTrim == 0){
	print SH "cp \$fastq/*.fastq.gz \$workdir/\$project.fastq/\n";
    }
    print SH "cd \$workdir\n";
}
if ($table eq ""){
    print SH "\n";
    print SH "mkdir \$workdir/readBased\n";
    print SH "\# Pull out the individual reads for read-based analyses.\n";
    print SH "for f in $project.fastq/\*.f*q.gz\n";
    print SH "do\n";
    print SH "\t\# Determine the sample name based on the fastq file name.\n";
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
	    print SH "\tgrep \^\@ -A 1 \$project.fastq/\$sample\\_trimmed.fq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > readBased/\$sample.justReads.fastq\n";
	}elsif ($runTrim == 0){
	    print SH "\tgrep \^\@ -A 1 \$project.fastq/\$sample.fastq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > readBased/\$sample.justReads.fastq\n";
	}
    } elsif ($barcodes ne "") {
	print SH "\tgrep \^\@ -A 1 \$project.fastq/\$sample.fastq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > readBased/\$sample.justReads.fastq\n";
    }
    print SH "\n\t# Combine duplicate reads and add read counts to each line.\n";
    print SH "\tsort readBased\/\$sample.justReads.fastq | uniq -c | sort -nr > readBased\/\$sample.justReads.uniqCounts.txt\n";
    print SH "\tgzip \$project.fastq/\$sample*q &\n";
    print SH "\tgzip readBased/\$sample.justReads.fastq &\n";
    print SH "done\n";
    print SH "\n";
    print SH "cd \$workdir/readBased\n";
    if ($barcodes ne ""){
	$type = "notUMId.";
	print SH "\# Process UMI sequences to split out reads into two sets, UMId and notUMId.\n";
	print SH "for f in \*.justReads.uniqCounts.txt\n";
	print SH "do\n";
	print SH "\t\# Extract the sample name from the file name.\n";
	print SH "\tfile=\$\(basename \$f\)\n";
	print SH "\techo \$file\n";
	print SH "\tsample=\$\{file\%\%.justReads.uniqCounts.txt\}\n";
	print SH "\techo \$sample\n";
	print SH "\n\t\# Generate two files, and remove UMI sequences from reads.\n";
	print SH "\t# notUMId: A six nucleotide Unique Molecular Identifier (UMI) is removed from each read, and the numbers of reads associated with each UMI are summed to generate the raw read count for a particular sample.  For the UMId files, the number of unique UMIs is counted, regardless of how frequently each UMI is seen.\n";
	print SH "\t\(perl \$pipeline\/processUMIs.pl \$sample.justReads.uniqCounts.txt\) \>\& \$sample.justReads.uniqCounts.UMI.dist.error.log\n"; 
	print SH "done\n";
    }
    print SH "\n# Remove adapter sequences from raw read counts.\n";
    print SH "for f in \*.justReads.uniqCounts*txt\n";
    print SH "do\n";
    print SH "\t# Adapter reads are identified as containing the substring GTCCGACGATC followed by 3-5 random nucleotides and removed from analysis.\n";
    print SH "\tperl \$pipeline\/deleteAdapterReads.pl \$f > \$f.noAdapter.txt\n";
    print SH "done\n";
    print SH "\n# Make tables from the reads.\n";
    print SH "# Counts are tabulated for each read and sample. rawCounts has raw un-normalized data.\n";
    #    print SH "perl \$pipeline/combineCounts.pl justReads.uniqCounts.".$type."txt.noAdapter.txt > \$project.noAdapter.".$type."normCounts.table.txt & \n";
    print SH "perl \$pipeline/combineCountsNotNorm.pl justReads.uniqCounts.".$type."txt.noAdapter.txt 1 > \$project.noAdapter.".$type."rawCounts.table.txt & \n";
    $table = "$project.noAdapter.".$type."rawCounts.table.txt";
    print SH "wait\n";
}
print SH "# Using Table of raw counts $table for the rest of the pipeline.\n";

my $tablePrefix = "";
if (basename($table) =~ /^([\w\_\.\-]+).txt$/){
    $tablePrefix = $1;
}
print SH "\n# Normalize raw read counts to counts per million, removing reads with fewer than $minSum counts across all samples.\n";
print SH "perl \$pipeline/normalizeRawCountsToCPM.pl $table $minSum > $tablePrefix.normCounts.txt\n";

print SH "# Make a fasta formatted file from the reads in the normCounts table\n";
print SH "awk \' \$1 \!\~ \"Read\" \{printf \"\>\%s.\%d\\n\%s\\n\",\$1,\$NF,\$1\}\' $tablePrefix.normCounts.txt > $tablePrefix.reads.fa\n";
print SH "# Reads were blasted (blastn-short) against custom blast databases created from all processed miRNAs and from the most recent RNA world databases, using the sequences appropriate for the organism (human or mouse). \n";
if ($organism eq "human"){
    print SH "# Blast reads against human miRNAs\n";
    print SH "blastn -task blastn-short -query $tablePrefix.reads.fa -db $questSporos/usefulFiles/humanMiRNAs.nr -outfmt 7 -num_threads $threads > $tablePrefix.miRNAs.nr.$organism.blast.txt \& \n";
    print SH "\# Blast reads against human RNA world.\n";
    print SH "blastn -task blastn-short -query $tablePrefix.reads.fa -db $questSporos/usefulFiles/human_and_virus_vbrc_all_cali_annoDB_sep2015.proc -outfmt 7 -num_threads $threads > $tablePrefix.RNAworld.$organism.blast.txt &\n";
} elsif ($organism eq "mouse"){
    print SH "\# Blast reads against mouse miRNAs\n";
    print SH "blastn -task blastn-short -query $tablePrefix.reads.fa -db $questSporos/usefulFiles/mouseMiRNAs.nr -outfmt 7 -num_threads $threads > $tablePrefix.miRNAs.nr.$organism.blast.txt \& \n";
    print SH "\# Blast reads against mouse RNA world.\n";
    print SH "blastn -task blastn-short -query $tablePrefix.reads.fa -db $questSporos/usefulFiles/mouse_annoDB_piRNA_all_oct2013.proc -outfmt 7 -num_threads $threads > $tablePrefix.RNAworld.$organism.blast.txt &\n";
}
print SH "wait\n";
print SH "\n\# After blast runs are complete, resume analysis.\n";
print SH "\n\# Analyze blast data.\n";
print SH "# Filter blast results\n";
print SH "# Processed miRNA hits were filtered for 100% identity, and a match length of at least 18 bp.\n";
print SH "awk \' \$3 == 100 && \$8 > \$7 && \$10 > \$9 && \$4 >= 18 && \$9 < 10 \' $tablePrefix.miRNAs.nr.$organism.blast.txt  > $tablePrefix.miRNAs.nr.blast.filtered.txt\n";
print SH "# RNAworld hits were filtered for at least 95% identity, where the match is within the first 9 bp of the sequencing reads.\n";
print SH "awk \' \$3 >= 95 && \$8 > \$7 && \$10 > \$9 && \$7 < 10\'  $tablePrefix.RNAworld.$organism.blast.txt  > $tablePrefix.RNAworld.blast.filtered.txt\n";
print SH "\n # Add seeds, toxicities, blast hits to normalized counts table.\n";
print SH "perl $questSporos/usefulFiles/addAllToxicities.pl $questSporos/usefulFiles/6mer\\ Seed\\ toxes\\ new.txt $tablePrefix.normCounts.txt > $tablePrefix.normCounts.withTox.txt\n";
print SH "perl \$pipeline/addBLASTresults.pl $tablePrefix.normCounts.withTox.txt $tablePrefix.miRNAs.nr.blast.filtered.txt > $tablePrefix.normCounts.withTox.withMiRNAblast.txt\n";
print SH "perl \$pipeline/addBLASTresults.pl $tablePrefix.normCounts.withTox.withMiRNAblast.txt $tablePrefix.RNAworld.blast.filtered.txt > $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt\n";
print SH "perl \$pipeline/truncateBLASTresults.pl $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt > $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt\n";

print SH "\n# Pull out smaller subsets of these tables for Excel.\n";
#print SH "perl \$pipeline/thresholdNormTotal.pl $tablePrefix.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.txt 20 > $tablePrefix.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.minSum20.txt\n";
print SH "perl \$pipeline/thresholdNormTotal.pl $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt 6 > $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.minSum6.txt\n";
print SH "\n# This \"normCounts\" file is used as the basis for many further analyses in the Peter Lab.  For simplicity, we will rename it normCounts.\$project.txt \n";
print SH "ln -s $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt normCounts.\$project.txt\n";
print SH "\n# Collapse counts for seed sequences from the same RNA species.\n";
print SH "perl \$pipeline/collapseSpecies.pl normCounts.\$project.txt > collapsed.\$project.txt\n";
print SH "\n# Collapse counts for seed sequences, regardless of RNA species of origin.\n";
print SH "perl \$pipeline/collapseToxicityBins.pl collapsed.\$project.txt $organism $sequenceType\n";
print SH "\n# Expand seed counts for each samples (and the average of any replicates) into lines in a text file, for weblogo.\n";
print SH "perl \$pipeline/expandSequencesFromSeedKeyed.pl seedKeyed.$sequenceType.$organism.\$project.txt $sequenceType \$project\n";
print SH "\n# Expand tox counts for each sample (and the average of any replicates) into lines in a text file, for boxplots.\n";
print SH "perl \$pipeline/expandToxesFromSeedCollapsed.pl collapsed.\$project.txt $organism $sequenceType \$project\n";
print SH "\n# Run weblogo.\n";
print SH "module load python/anaconda3.6\n";
print SH "source activate /projects/p20742/envs/weblogo-py38\n";
print SH "for f in seedAnalysis*.txt\ndo\n";
print SH "\techo \$f\n";
print SH "\tsample=\$\{f\%\%.txt\}\n";
print SH "\tweblogo -f \$f -A rna -U probability -F pdf -o \$sample.pdf --color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t--size large --ylabel Probability --logo-font Helvetica-Extra-Bold --number-interval 1\n";
print SH "done \n";
print SH "source deactivate\n";


if ($comparisons ne ""){
    print SH "module load R/3.2.1\n";
    print SH "\n# Run EdgeR on the comparisons file to identify differentially abundant reads.\n";
    print SH "# (Note that assembly is not actually used here, because these are reads, not gene IDs)\n";
    print SH "perl -pe \"s\/.justReads.uniqCounts.".$type."txt\/\/g\" $table | perl -pe \"s/.noAdapter.txt//g\" | perl -pe \"s\/\_R1\_001\/\/g\" | perl -pe \"s/.geo.txt//g\" > $tablePrefix.relabeled.txt\n";
    print SH "Rscript \$pipeline/runEdgeRrnaSeq.reads.R --countFile=$tablePrefix.relabeled.txt --numCores=$threads --runMDS=1 --filterOff=0 --assembly=hg38.mp\n";
    print SH "Rscript \$pipeline/runEdgeRrnaSeq.reads.R --countFile=$tablePrefix.relabeled.txt --numCores=$threads --runMDS=0 --filterOff=0 --assembly=hg38.mp --comparisonFile=$comparisons\n";
    print SH "\n# Annotate each EdgeR result file with toxicities and blast results.\n";
    print SH "for edgeR in *txt.edgeR.txt\n";
    print SH "do\n";
    print SH "\techo \$edgeR\n";
    print SH "\tcomp=\$\{edgeR\%\%.$tablePrefix.relabeled.txt.edgeR.txt\}\n";
    print SH "\techo \$comp\n";
    print SH "\tperl -pe \"s\/\\\"\/\/g\" \$edgeR > \$comp.edgeR.txt\n";
    print SH "\n\t# Add seed toxicity and source sequences to edgeR results.\n";
    print SH "\tperl $questSporos/usefulFiles/addAllToxicities.pl $questSporos/usefulFiles/6mer\\ Seed\\ toxes\\ new.txt \$comp.edgeR.txt > \$comp.edgeR.withTox.txt\n";
    print SH "\tperl \$pipeline/addBLASTresults.pl \$comp.edgeR.withTox.txt $tablePrefix.miRNAs.nr.blast.filtered.txt > \$comp.edgeR.withTox.withMiRNAblast.txt\n";
    print SH "\tperl \$pipeline/addBLASTresults.pl \$comp.edgeR.withTox.withMiRNAblast.txt $tablePrefix.RNAworld.blast.filtered.txt > \$comp.edgeR.withTox.withMiRNAandRNAworld.blast.txt\n";
    print SH "\tperl \$pipeline/truncateBLASTresults.pl \$comp.edgeR.withTox.withMiRNAandRNAworld.blast.txt > \$comp.edgeR.withTox.withMiRNAandRNAworld.blast.trunc.txt\n";
    print SH "\n\t# Rename file for simplicity.\n";
    print SH "\tln -s \$comp.edgeR.withTox.withMiRNAandRNAworld.blast.trunc.txt normCounts.\$comp.edgeR.txt\n";
    print SH "\n\t# Separate out reads that are either up or down regulated at an adjp <= 0.05 and a logFC > 0.585\n";
    print SH "\tperl \$pipeline/thresholdDifferential.pl normCounts.\$comp.edgeR.txt adjp 0.585 0.05\n";
    print SH "\n\tfor de in up dn\n";
    print SH "\tdo\n";
    print SH "\techo \$de\n";
    print SH "\n\t# Collapse counts for seed sequences from the same RNA species.\n";
    print SH "\tperl \$pipeline/collapseSpecies.pl normCounts.\$comp.edgeR.\$de.txt > collapsed.\$comp.edgeR.\$de.txt\n";
    print SH "\n\t# Collapse counts for seed sequences, regardless of RNA species of origin.\n";
    print SH "\tperl \$pipeline/collapseToxicityBins.pl collapsed.\$comp.edgeR.\$de.txt $organism $sequenceType\n";
    print SH "\n\t# Expand seed counts for each sample (and the average of any replicates) into lines in a text file, for weblogo.\n";
    print SH "\tperl \$pipeline/expandSequencesFromSeedKeyed.pl seedKeyed.$sequenceType.$organism.\$comp.edgeR.\$de.txt $sequenceType \$comp.\$de\n";
    print SH "\n\t# Expand tox counts for each sample (and the average of any replicates) into lines in a text file, for boxplots.\n";
    print SH "\tperl \$pipeline/expandToxesFromSeedCollapsed.pl collapsed.\$comp.edgeR.\$de.txt $organism $sequenceType \$comp.\$de\n";
    print SH "\tdone\n";
    print SH "\n\t# Group together seedAnalysis sequences from samples from the same comparison grouping.\n";
    print SH "\tperl \$pipeline/prepareComparisonWeblogoSequences.pl $comparisons \$project\n";
    print SH "\n\t# Run weblogo.\n";
    print SH "\tmodule load python/anaconda3.6\n";
    print SH "\tsource activate /projects/p20742/envs/weblogo-py38\n";
    print SH "\tfor f in seedAnalysis*\$comp*txt\n\tdo\n";
    print SH "\techo \$f\n";
    print SH "\tsample=\$\{f\%\%.txt\}\n";
    print SH "\tsort \$f > \$sample.sorted.txt\n";
    print SH "\tweblogo -f \$sample.sorted.txt -A rna -U probability -F pdf -o \$sample.pdf --color \"\#CC0000\" G guanine --color \"\#FFB302\" C cytosine --color \"\#0100CC\" A adenine --color \"\#01CC00\" U uracil \\\n\t --size large --ylabel Probability --logo-font Helvetica-Extra-Bold --number-interval 1\n";
    print SH "\tdone \n";
    print SH "\tsource deactivate\n";
    print SH "\tdone \n";
    print SH "done\n";
}
