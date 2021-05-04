#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -m a
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1069/GSE63501.autoRepeat
#SBATCH -o "%x.o%j"
#SBATCH --job-name=GSE63501
#SBATCH --nodes=1
#SBATCH -n 24

pipeline=/projects/b1069//smallRNAscripts/
project=GSE63501
workdir=/projects/b1069/GSE63501.autoRepeat
fastq=/projects/b1069/GSE63501.autoRepeat/fastq
currentDate=$(date +%F)

module load python/anaconda
module load R/3.2.1
module load bedtools
module load blast/2.7.1

cd $workdir
mkdir $workdir/$project.fastq


# Trim the reads to remove any standard illumina adaptors.
for f in $fastq/*.fastq.gz
do
	echo $f
	perl /projects/p20742/tools/bin//trim_galore $f --length 6 &
done

# Wait for trimming to finish.
wait
# Move trimmed files into $project.fastq directory.
mv *trimmed.fq.gz $workdir/$project.fastq/

mkdir $workdir/readBased
# Pull out the individual reads for read-based analyses.
for f in GSE63501.fastq/*.f*q.gz
do
	# Determine the sample name based on the fastq file name.
	file=$(basename $f)
	echo $file
	sample=${file%%.gz}
	sample=${sample%%_trimmed.fq}
	sample=${sample%%.fastq}

	# Print the sample name for the fastq file.
	echo $sample
	gunzip $f

	# Pull out just the read sequences, leaving quality scores behind.
	grep ^@ -A 1 $project.fastq/$sample\_trimmed.fq | grep -e ^A -e ^T -e ^G -e ^C -e ^N > readBased/$sample.justReads.fastq

	# Combine duplicate reads and add read counts to each line.
	sort readBased/$sample.justReads.fastq | uniq -c | sort -nr > readBased/$sample.justReads.uniqCounts.txt
	gzip $project.fastq/$sample*q &
	gzip readBased/$sample.justReads.fastq &
done

cd $workdir/readBased

# Make tables from the reads.
# Counts are tabulated for each read and sample and two files are generated.  normCounts has the counts normalized per million reads in the column (sample) sum.  rawCounts has raw un-normalized data.
perl $pipeline/combineCounts.pl justReads.uniqCounts.txt > $project.allreads.normCounts.table.txt & 
perl $pipeline/combineCountsNotNorm.pl justReads.uniqCounts.txt > $project.allreads.rawCounts.table.txt & 
wait
# Adapter reads are identified as containing the substring GTCCGACGATC followed by 3-5 random nucleotides and removed from analysis.
perl $pipeline/deleteAdapterReads.pl $project.allreads.rawCounts.table.txt > $project.noAdapter.rawCounts.table.txt
perl $pipeline/deleteAdapterReads.pl $project.allreads.normCounts.table.txt > $project.noAdapter.normCounts.table.txt
# Make a fasta formatted file from the reads in the normCounts table
awk ' $1 !~ "Read" {printf ">%s.%d\n%s\n",$1,$NF,$1}' $project.noAdapter.normCounts.table.txt > $project.noAdapter.fa
# Reads were blasted (blastn-short) against custom blast databases created from all processed miRNAs and from the most recent RNA world databases, using the sequences appropriate for the organism (human or mouse). 
# Blast reads against human miRNAs
blastn -task blastn-short -query $project.noAdapter.fa -db /projects/b1069//usefulFiles/humanMiRNAs.nr -outfmt 7 -num_threads 24 > $project.noAdapter.miRNAs.nr.human.blast.txt & 
# Blast reads against human RNA world.
blastn -task blastn-short -query $project.noAdapter.fa -db /projects/b1069//usefulFiles/human_and_virus_vbrc_all_cali_annoDB_sep2015.proc -outfmt 7 -num_threads 24 > $project.noAdapter.RNAworld.human.blast.txt &
wait

# After blast runs are complete, resume analysis.

# Analyze blast data.
# Filter blast results
# Processed miRNA hits were filtered for 100% identity, and a match length of at least 18 bp.
awk ' $3 == 100 && $8 > $7 && $10 > $9 && $4 >= 18 && $9 < 10 ' $project.noAdapter.miRNAs.nr.human.blast.txt  > $project.noAdapter.miRNAs.nr.blast.filtered.txt
# RNAworld hits were filtered for at least 95% identity, where the match is within the first 9 bp of the sequencing reads.
awk ' $3 >= 95 && $8 > $7 && $10 > $9 && $7 < 10'  $project.noAdapter.RNAworld.human.blast.txt  > $project.noAdapter.RNAworld.blast.filtered.txt
for c in normCounts rawCounts
do
	echo $c
	perl /projects/b1069//usefulFiles/addAllToxicities.pl /projects/b1069//usefulFiles/6mer\ Seed\ toxes\ new.txt $project.noAdapter.$c.table.txt > $project.noAdapter.$c.table.withTox.txt
	perl $pipeline/addBLASTresults.pl $project.noAdapter.$c.table.withTox.txt $project.noAdapter.miRNAs.nr.blast.filtered.txt > $project.noAdapter.$c.table.withTox.withMiRNAblast.txt
	perl $pipeline/addBLASTresults.pl $project.noAdapter.$c.table.withTox.withMiRNAblast.txt $project.noAdapter.RNAworld.blast.filtered.txt > $project.noAdapter.$c.table.withTox.withMiRNAandRNAworld.blast.txt
	perl $pipeline/truncateBLASTresults.pl $project.noAdapter.$c.table.withTox.withMiRNAandRNAworld.blast.txt > $project.noAdapter.$c.table.withTox.withMiRNAandRNAworld.blast.trunc.txt
done

# Pull out smaller subsets of these tables for Excel.
perl $pipeline/thresholdNormTotal.pl $project.noAdapter.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.txt 20 > $project.noAdapter.rawCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.minSum20.txt
perl $pipeline/thresholdNormTotal.pl $project.noAdapter.normCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.txt 6 > $project.noAdapter.normCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.minSum6.txt

# Collapse counts for seed sequences from the same RNA species.
perl $pipeline/collapseSpecies.pl $project.noAdapter.normCounts.table.withTox.withMiRNAandRNAworld.blast.trunc.minSum6.txt > $project.normCounts.seedCollapsed.txt

# Collapse counts for seed sequences, regardless of RNA species of origin.
perl $pipeline/collapseToxicityBins.pl $project.normCounts.seedCollapsed.txt human

# Expand seed counts for each samples (and the average of any replicates) into lines in a text file, for weblogo.
perl $pipeline/expandSequencesFromSeedKeyed.pl $project.normCounts.seedCollapsed.human.seedKeyed.txt

# Run weblogo.
module load python/anaconda3.6
source activate /projects/p20742/envs/weblogo-py38
for f in *seeds.1000.txt
do
echo $f
	sample=${f%%.txt}
	weblogo -f $f -A rna -U probability -F pdf -o $sample.pdf --color "#CC0000" G guanine --color "#FFB302" C cytosine --color "#0100CC" A adenine --color "#01CC00" U uracil --size large --ylabel Probability --logo-font Helvetica-Extra-Bold --number-interval 1
done 
source deactivate
