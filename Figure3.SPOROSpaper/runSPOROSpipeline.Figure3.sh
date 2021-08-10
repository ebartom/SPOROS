#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -m a
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1069/Figure3.SPOROSpaper/
#SBATCH -o "%x.o%j"
#SBATCH --job-name=Figure3.SPOROS
#SBATCH --nodes=1
#SBATCH -n 24

pipeline=/projects/b1069//smallRNAscripts/
project=Figure3
workdir=/projects/b1069/Figure3.SPOROSpaper/
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
cd $workdir/prelimAnalysis
cp /projects/b1069/Figure3.SPOROSpaper/input/GSE63501.publishedRawCounts.minSum0.txt $workdir/prelimAnalysis

# Using table of raw counts /projects/b1069/Figure3.SPOROSpaper/input/GSE63501.publishedRawCounts.minSum0.txt for the rest of the pipeline.
cd $workdir/prelimAnalysis
tablePrefix=GSE63501.publishedRawCounts.minSum0
# Remove reads / rows with a count sum < the number of samples in the table. The last column in the table should be the row total.
awk  ' $NF >= (NF-2) ' /projects/b1069/Figure3.SPOROSpaper/input/GSE63501.publishedRawCounts.minSum0.txt > $tablePrefix.minSumN.txt

# Make a fasta formatted file from the reads in the table
awk ' $1 !~ "Read" {printf ">%s.%d\n%s\n",$1,$NF,$1}' GSE63501.publishedRawCounts.minSum0.minSumN.txt > $tablePrefix.reads.fa

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
perl -pe "s/GSE63501.publishedRawCounts.minSum0.RNAworld.blast.filtered.txt/RNAworld/g"  $workdir/prelimAnalysis/$tablePrefix.normCounts.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe "s/GSE63501.publishedRawCounts.minSum0.miRNAs.nr.blast.filtered.txt/miRNA/g" > $workdir/totalCounts/A_normCounts.$project.txt
# Also copy over rawCounts version.
perl -pe "s/GSE63501.publishedRawCounts.minSum0.RNAworld.blast.filtered.txt/RNAworld/g"  $workdir/prelimAnalysis/$tablePrefix.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe "s/GSE63501.publishedRawCounts.minSum0.miRNAs.nr.blast.filtered.txt/miRNA/g" > $workdir/totalCounts/A_rawCounts.$project.txt
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
			sample=${sample%%.human.Figure3.txt}
			echo $sample | cat > test.$sample.tox
			cat $d >> test.$sample.tox
		done
	# If there are no averaged files, just use all of the individual ones.
	else
		for d in D_toxAnalysis.*.txt
		do
			echo $d
			sample=${d##D_toxAnalysis.}
			sample=${sample%%.human.Figure3.txt}
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

# Running differential read analysis.

mkdir $workdir/differential
module load R/3.2.1
cd $workdir/prelimAnalysis

# Run EdgeR on the comparisons file to identify differentially abundant reads.
# (Note that assembly is not actually used here, because these are reads, not gene IDs)
perl -pe "s/.justReads.uniqCounts.txt//g" /projects/b1069/Figure3.SPOROSpaper/input/GSE63501.publishedRawCounts.minSum0.txt | perl -pe "s/.noAdapter.txt//g" | perl -pe "s/_S\d+_R1_001//g" | perl -pe "s/.geo.txt//g" > $tablePrefix.relabeled.txt
Rscript $pipeline/runEdgeRrnaSeq.reads.R --countFile=$tablePrefix.relabeled.txt --numCores=24 --runMDS=1 --filterOff=0 --assembly=hg38.mp
Rscript $pipeline/runEdgeRrnaSeq.reads.R --countFile=$tablePrefix.relabeled.txt --numCores=24 --runMDS=0 --filterOff=0 --assembly=hg38.mp --comparisonFile=/projects/b1069/Figure3.SPOROSpaper/input/comparisons.csv

# Annotate each EdgeR result file with toxicities and blast results.
for edgeR in *relabeled.txt.edgeR.txt
do
	echo $edgeR
	comp=${edgeR%%.GSE63501.publishedRawCounts.minSum0.relabeled.txt.edgeR.txt}
	echo $comp
	cd $workdir/prelimAnalysis/
	perl -pe "s/\"//g" $edgeR > $comp.txt

	# Add seed toxicity and source sequences to edgeR results.
	perl /projects/b1069//usefulFiles/addAllToxicities.pl /projects/b1069//usefulFiles/speciesToxes.txt $comp.txt > $comp.withTox.txt
	perl $pipeline/addBLASTresults.pl $comp.withTox.txt $tablePrefix.miRNAs.nr.blast.filtered.txt > $comp.withTox.withMiRNAblast.txt
	perl $pipeline/addBLASTresults.pl $comp.withTox.withMiRNAblast.txt $tablePrefix.RNAworld.blast.filtered.txt > $comp.withTox.withMiRNAandRNAworld.blast.txt
	perl $pipeline/truncateBLASTresults.pl $comp.withTox.withMiRNAandRNAworld.blast.txt > $comp.withTox.withMiRNAandRNAworld.blast.trunc.txt

	# Rename file for simplicity and get organized.
	mkdir $workdir/differential/$comp
	perl -pe "s/GSE63501.publishedRawCounts.minSum0.RNAworld.blast.filtered.txt/RNAworld/g"  $workdir/prelimAnalysis/$comp.withTox.withMiRNAandRNAworld.blast.trunc.txt | perl -pe "s/GSE63501.publishedRawCounts.minSum0.miRNAs.nr.blast.filtered.txt/miRNA/g" > $workdir/differential/$comp/A_normCounts.$comp.all.txt
	cd $workdir/differential/$comp
	# Create a subdirectory for both an adjusted p-value comparison and a p-value comparison, and then run both of them within the relevant sub-directory
	for stat in adjp pvalue
	do
		mkdir $workdir/differential/$comp/$stat

		# Separate out reads that are differentially regulated at a $stat <= 0.05 and a absolute logFC > 0.585
		cd $workdir/differential/$comp/$stat
		perl $pipeline/thresholdDifferential2.pl ../A_normCounts.$comp.all.txt $stat 0.585 0.05
		# This generates files A_normCounts.$comp.diff.txt with all differentially regulated reads.


		# Collapse counts for seed sequences from the same RNA species.
		perl $pipeline/collapseSpecies.pl A_normCounts.$comp.diff.txt > B_collapsed.$comp.diff.txt

		# Guided by comparisons file, find the average counts across the two groups indicated for the comparison $comp, for differential reads
		perl $pipeline/averageWithinGroups.pl /projects/b1069/Figure3.SPOROSpaper/input/comparisons.csv B_collapsed.$comp.diff.txt $comp > B_collapsed.$comp.diff.avgd.txt
		# Split downstream analysis depending on whether the reads from miRNAs or sRNAs.
		for seqtype in sRNA miRNA
		do
			mkdir $workdir/differential/$comp/$stat/$seqtype
			cd $workdir/differential/$comp/$stat/$seqtype

			# Collapse counts for differential and up-regulated read sequences, regardless of RNA species of origin, as long as they are $seqtype.
			perl $pipeline/collapseToxicityBins.pl ../B_collapsed.$comp.diff.avgd.txt human $seqtype
			# This generates files C_binned.$seqtype.human.$comp.diff.avgd.txt and intermediary file Int_seedKeyed.$seqtype.human.$comp.diff.avgd.txt.

			# Subtract read counts in Group1 from Group2 to find the difference between the two.
			perl $pipeline/deltaBetweenGroups.pl Int_seedKeyed.$seqtype.human.$comp.diff.avgd.txt
			# This generates files Int_seedKeyed.$seqtype.human.$comp.delta.up.txt and Int_seedKeyed.$seqtype.human.$comp.delta.up.txt 

			# Expand delta seed counts for seeds that went either up or down into lines in a text file, for weblogo.
			perl $pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.$seqtype.human.$comp.delta.up.txt $seqtype up
			perl $pipeline/expandSequencesFromSeedKeyed.pl Int_seedKeyed.$seqtype.human.$comp.delta.dn.txt $seqtype dn

			# Expand tox counts for each delta up or delta down file into lines in a text file, for boxplots.
			perl $pipeline/expandToxesFromSeedKeyed.pl Int_seedKeyed.$seqtype.human.$comp.delta.up.txt up
			perl $pipeline/expandToxesFromSeedKeyed.pl Int_seedKeyed.$seqtype.human.$comp.delta.dn.txt dn

			# Run weblogo.
			module load python/anaconda3.6
			source activate /projects/p20742/envs/weblogo-py38
			for f in E_seedAnalysis*Delta*txt
			do
				echo $f
				sample=${f%%.txt}
				weblogo -f $sample.txt -A rna --show-xaxis NO --show-yaxis NO -U probability -F png_print -o $sample.png \
					--color "#CC0000" G guanine --color "#FFB302" C cytosine --color "#0100CC" A adenine --color "#01CC00" U uracil \
					--size large --logo-font Helvetica-Extra-Bold --ylabel Probability --number-interval 1
				done
			source deactivate

			# Compile the delta toxAnalysis files into one multi-column summary.
			for d in D_toxAnalysis.delta.*.txt
			do
				echo $d
				sample=${d##D_toxAnalysis.}
				sample=${sample%%.human.Figure3.txt}
				echo $sample | cat > test.$sample.tox
				cat $d >> test.$sample.tox
			done
			paste test.*.tox > D_toxAnalysis.combined.$seqtype.txt
			rm test.*.tox

			# Plot a very basic box plot summarizing the tox distributions in the different samples.
			Rscript /projects/b1069/smallRNAscripts/plotToxBoxPlot.R --toxFile=D_toxAnalysis.combined.$seqtype.txt
			done
		done
	done
done
