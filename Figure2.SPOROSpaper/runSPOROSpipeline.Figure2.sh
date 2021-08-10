#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -m a
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1069/Figure2.SPOROSpaper
#SBATCH -o "%x.o%j"
#SBATCH --job-name=Figure2.SPOROS
#SBATCH --nodes=1
#SBATCH -n 24

pipeline=/projects/b1069//smallRNAscripts/
project=Figure2
workdir=/projects/b1069/Figure2.SPOROSpaper
fastq=/projects/b1069/Figure2.SPOROSpaper/input/fastq/
currentDate=$(date +%F)

module load perl/5.16
module load python/anaconda
module load R/3.2.1
module load bedtools/2.29.2
module load blast/2.7.1

cd $workdir
# Starting analysis for SPOROS pipeline.
date

mkdir $workdir/prelimAnalysis
mkdir $workdir/$project.fastq


# Barcode files for demultiplexing found.  /projects/b1069/Figure2.SPOROSpaper/input/HCT116_Drosha_Dicer_rep1.barcodes.txt /projects/b1069/Figure2.SPOROSpaper/input/HCT116_Drosha_Dicer_rep2.barcodes.txt

# Each barcode needs to match a fastq file in the /projects/b1069/Figure2.SPOROSpaper/input/fastq/ directory, with the same file prefix.
mkdir $workdir/demultiplex
cd $workdir/demultiplex

cp /projects/b1069/Figure2.SPOROSpaper/input/HCT116_Drosha_Dicer_rep1.barcodes.txt $workdir/demultiplex/ 
# Demultiplexing samples for sample group HCT116_Drosha_Dicer_rep1
perl /projects/p20742/tools/bin//trim_galore $fastq/HCT116_Drosha_Dicer_rep1.fastq.gz --length 12 --dont_gzip -o $workdir/demultiplex/ 
perl $pipeline/splitFastQwithTable.pl HCT116_Drosha_Dicer_rep1_trimmed.fq HCT116_Drosha_Dicer_rep1.barcodes.txt 0
mv *.fastq $workdir/$project.fastq/
gzip $workdir/$project.fastq/*.fastq &
gzip $workdir/demultiplex/HCT116_Drosha_Dicer_rep1*.fq &

cp /projects/b1069/Figure2.SPOROSpaper/input/HCT116_Drosha_Dicer_rep2.barcodes.txt $workdir/demultiplex/ 
# Demultiplexing samples for sample group HCT116_Drosha_Dicer_rep2
perl /projects/p20742/tools/bin//trim_galore $fastq/HCT116_Drosha_Dicer_rep2.fastq.gz --length 12 --dont_gzip -o $workdir/demultiplex/ 
perl $pipeline/splitFastQwithTable.pl HCT116_Drosha_Dicer_rep2_trimmed.fq HCT116_Drosha_Dicer_rep2.barcodes.txt 0
mv *.fastq $workdir/$project.fastq/
gzip $workdir/$project.fastq/*.fastq &
gzip $workdir/demultiplex/HCT116_Drosha_Dicer_rep2*.fq &
wait

# Finished de-multiplexing.
date
cd $workdir

# Building counts table from raw reads.
mkdir $workdir/prelimAnalysis
cd $workdir/prelimAnalysis
# Pull out the individual reads for read-based analyses.
for f in $workdir/Figure2.fastq/*.f*q.gz
do
	# Determine the sample name based on the fastq file name.
	echo $f
	file=$(basename $f)
	echo $file
	sample=${file%%.gz}
	sample=${sample%%_trimmed.fq}
	sample=${sample%%.fastq}

	# Print the sample name for the fastq file.
	echo $sample
	gunzip $f

	# Pull out just the read sequences, leaving quality scores behind.
	grep ^@ -A 1 $workdir/$project.fastq/$sample.fastq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > $workdir/prelimAnalysis/$sample.justReads.fastq

	# Combine duplicate reads and add read counts to each line.
	sort $workdir/prelimAnalysis/$sample.justReads.fastq | uniq -c | sort -nr > $workdir/prelimAnalysis/$sample.justReads.uniqCounts.txt
	gzip $workdir/$project.fastq/$sample*q &
	gzip $workdir/prelimAnalysis/$sample.justReads.fastq &
done

cd $workdir/prelimAnalysis

# We assume that samples with barcodes also have UMI sequences of the pattern 5'-NNNNreadNN-3'.  If this is not the case, fastq files should be demultiplexed before starting SPOROS, and it should be run without demultiplexing barcodes.

# Process UMI sequences to split out reads into two sets, UMId and notUMId.
for f in *.justReads.uniqCounts.txt
do
	# Extract the sample name from the file name.
	file=$(basename $f)
	echo $file
	sample=${file%%.justReads.uniqCounts.txt}
	echo $sample

	# Generate two files, and remove UMI sequences from reads.
	# notUMId: A six nucleotide Unique Molecular Identifier (UMI) is removed from each read, and the numbers of reads associated with each UMI are summed to generate the raw read count for a particular sample.
	# For the UMId files, the number of unique UMIs is counted, regardless of how frequently each UMI is seen.
	perl $pipeline/processUMIs.pl $sample.justReads.uniqCounts.txt
done

# Remove adapter sequences from raw read counts.
for f in *.justReads.uniqCounts*txt
do
	# Adapter reads are identified as containing the substring TCCGACGATC followed by 3-5 random nucleotides and removed from analysis.
	perl $pipeline/deleteAdapterReads.pl $f > $f.noAdapter.txt
done

# Make tables from the reads.
# Counts are tabulated for each read and sample. The rawCounts table has raw un-normalized data.
perl $pipeline/combineCountsNotNorm.pl justReads.uniqCounts.notUMId.txt.noAdapter.txt 1 > $project.noAdapter.notUMId.rawCounts.table.txt 

# Using table of raw counts Figure2.noAdapter.notUMId.rawCounts.table.txt for the rest of the pipeline.
cd $workdir/prelimAnalysis
tablePrefix=Figure2.noAdapter.notUMId.rawCounts.table
# Remove reads / rows with a count sum < the number of samples in the table. The last column in the table should be the row total.
awk  ' $NF >= (NF-2) ' Figure2.noAdapter.notUMId.rawCounts.table.txt > $tablePrefix.minSumN.txt

# Make a fasta formatted file from the reads in the table
awk ' $1 !~ "Read" {printf ">%s.%d\n%s\n",$1,$NF,$1}' Figure2.noAdapter.notUMId.rawCounts.table.minSumN.txt > $tablePrefix.reads.fa

# Reads are blasted (blastn-short) against custom blast databases created from all processed miRNAs and from the most recent RNA world databases, using the sequences appropriate for the organism (human or mouse). 
# Blast reads against human miRNAs
blastn -task blastn-short -query $tablePrefix.reads.fa -db /projects/b1069//usefulFiles/humanMiRNAs.nr -outfmt 7 -num_threads 24 > $tablePrefix.miRNAs.nr.human.blast.txt & 
# Blast reads against human RNA world.
blastn -task blastn-short -query $tablePrefix.reads.fa -db /projects/b1069//usefulFiles/human_and_virus_vbrc_all_cali_annoDB_sep2015.proc -outfmt 7 -num_threads 24 > $tablePrefix.RNAworld.human.blast.txt &
wait

# After blast runs are complete, resume analysis.

# Analyze blast data.
# Filter blast results
# Processed miRNA hits were filtered for 100% identity, and a match length of at least 18 bp.
awk ' $3 == 100 && $8 > $7 && $10 > $9 && $4 >= 18 && $9 < 10 ' $tablePrefix.miRNAs.nr.human.blast.txt  > $tablePrefix.miRNAs.nr.blast.filtered.txt
# RNAworld hits were filtered for at least 95% identity, where the match is within the first 9 bp of the sequencing reads.
awk ' $3 >= 95 && $8 > $7 && $10 > $9 && $7 < 10'  $tablePrefix.RNAworld.human.blast.txt  > $tablePrefix.RNAworld.blast.filtered.txt

# Remove reads from table that have hits against RNAworld adapter sequences.
perl $pipeline/removeRNAworldAdapterSeq.pl $tablePrefix.minSumN.txt $tablePrefix.RNAworld.blast.filtered.txt > $tablePrefix.RNAworldClean.txt

# Normalize raw read counts to counts per million, removing reads with fewer than 2 counts across all samples. (Note that reads with fewer reads than the number of samples have already been removed.)
perl $pipeline/normalizeRawCountsToCPM.pl $tablePrefix.RNAworldClean.txt 2 > $tablePrefix.normCounts.txt

 # Add seeds, toxicities, blast hits to normalized counts table.
perl /projects/b1069//usefulFiles/addAllToxicities.pl /projects/b1069//usefulFiles/speciesToxes.txt $tablePrefix.normCounts.txt > $tablePrefix.normCounts.withTox.txt
perl $pipeline/addBLASTresults.pl $tablePrefix.normCounts.withTox.txt $tablePrefix.miRNAs.nr.blast.filtered.txt > $tablePrefix.normCounts.withTox.withMiRNAblast.txt
perl $pipeline/addBLASTresults.pl $tablePrefix.normCounts.withTox.withMiRNAblast.txt $tablePrefix.RNAworld.blast.filtered.txt > $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt
perl $pipeline/truncateBLASTresults.pl $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.txt > $tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt

 # Add seeds, toxicities, blast hits to raw counts table.
perl /projects/b1069//usefulFiles/addAllToxicities.pl /projects/b1069//usefulFiles/speciesToxes.txt $tablePrefix.RNAworldClean.txt > $tablePrefix.withTox.txt
perl $pipeline/addBLASTresults.pl $tablePrefix.withTox.txt $tablePrefix.miRNAs.nr.blast.filtered.txt > $tablePrefix.withTox.withMiRNAblast.txt
perl $pipeline/addBLASTresults.pl $tablePrefix.withTox.withMiRNAblast.txt $tablePrefix.RNAworld.blast.filtered.txt > $tablePrefix.withTox.withMiRNAandRNAworld.blast.txt
perl $pipeline/truncateBLASTresults.pl $tablePrefix.withTox.withMiRNAandRNAworld.blast.txt > $tablePrefix.withTox.withMiRNAandRNAworld.blast.trunc.txt

# This "normCounts" file is used as the basis for many further analyses in the Peter Lab, and is the first of the official output files for the SPOROS pipelien.  For simplicity, we will rename it A_normCounts.$project.txt and put it in a new directory with the other output files.
mkdir $workdir/totalCounts
perl -pe "s/Figure2.noAdapter.notUMId.rawCounts.table.RNAworld.blast.filtered.txt/RNAworld/g"  $workdir/prelimAnalysis/$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe "s/Figure2.noAdapter.notUMId.rawCounts.table.miRNAs.nr.blast.filtered.txt/miRNA/g" > $workdir/totalCounts/A_normCounts.$project.txt
# Also copy over rawCounts version.
perl -pe "s/Figure2.noAdapter.notUMId.rawCounts.table.RNAworld.blast.filtered.txt/RNAworld/g"  $workdir/prelimAnalysis/$tablePrefix.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe "s/Figure2.noAdapter.notUMId.rawCounts.table.miRNAs.nr.blast.filtered.txt/miRNA/g" > $workdir/totalCounts/A_rawCounts.$project.txt
cd $workdir/totalCounts

# Collapse counts for seed sequences from the same RNA species.
perl $pipeline/collapseSpecies.pl A_normCounts.$project.txt > B_collapsed.$project.txt

# Split the analysis by sequence type and generate a directory for each and run separately.
for seqtype in sRNA miRNA
do
	mkdir $workdir/totalCounts/$seqtype
	cd $workdir/totalCounts/$seqtype

	# Collapse counts for seed sequences, regardless of RNA species of origin, for all $seqtype.
	perl $pipeline/collapseToxicityBins.pl ../B_collapsed.$project.txt human $seqtype
	# This generates files C_binned.$seqtype.human.$project.txt and intermediary file Int_seedKeyed.$seqtype.human.$project.txt.

	# Expand tox counts for each sample (and the average of any replicates) into lines in a text file, for boxplots.
	perl $pipeline/expandToxesFromSeedCollapsed.pl ../B_collapsed.$project.txt human $seqtype $project
	# This generates files that start with the word D_toxAnalysis, one for each sample and average of replicate samples.

	# Expand seed counts for each samples (and the average of any replicates) into lines in a text file, for weblogo.
	perl $pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.$seqtype.human.$project.txt $seqtype $project
	# This generates files that start with the words E_seedAnalysis and F_seedExpand, one for each sample and average of replicate samples.

	# Run weblogo.
	module load python/anaconda3.6
	source activate /projects/p20742/envs/weblogo-py38
	for f in E_seedAnalysis*.txt
	do
		echo $f
		sample=${f%%.txt}
		weblogo -f $f -A rna --show-xaxis NO --show-yaxis NO -U probability -F png_print -o $sample.png \
			--color "#CC0000" G guanine --color "#FFB302" C cytosine --color "#0100CC" A adenine --color "#01CC00" U uracil \
			--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1
		done 
	source deactivate

	# Compile the toxAnalysis files into one multi-column summary.
	# If there are averaged files, use those.
	if ls | grep -q "avg"; then
		for d in D_toxAnalysis.*.avg.*.txt
		do
			echo $d
			sample=${d##D_toxAnalysis.}
			sample=${sample%%.human.Figure2.txt}
			echo $sample | cat > test.$sample.tox
			cat $d >> test.$sample.tox
		done
	# If there are no averaged files, just use all of the individual ones.
	else
		for d in D_toxAnalysis.*.txt
		do
			echo $d
			sample=${d##D_toxAnalysis.}
			sample=${sample%%.human.Figure2.txt}
			echo $sample | cat > test.$sample.tox
			cat $d >> test.$sample.tox
		done
	fi
	paste test.*.tox > D_toxAnalysis.combined.$seqtype.txt
	rm test.*.tox

	# Plot a very basic box plot summarizing the tox distributions in the different samples.
	Rscript /projects/b1069/smallRNAscripts/plotToxBoxPlot.R --toxFile=D_toxAnalysis.combined.$seqtype.txt

	# Analysis done for $seqtype.
	done
date

# Analysis done for total counts.
