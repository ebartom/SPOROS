#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 06:00:00
#SBATCH -m a
#SBATCH --mem=50000
#SBATCH --chdir=/projects/b1069/Figure3.SPOROSpaper/input/
#SBATCH -o "%x.o%j"
#SBATCH --job-name=Figure3.SPOROS.buildTable
#SBATCH --nodes=1
#SBATCH -n 12

pipeline=/projects/b1069//smallRNAscripts/
project=Figure3
workdir=/projects/b1069/Figure3.SPOROSpaper/input/
currentDate=$(date +%F)

# First download all data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63501
# Decompress tar file
 tar -xvf GSE63501_RAW.tar.gz

# Unzip data for each sample.
gunzip GSM* &
wait

# Sort each sample and concatenate.
for sample in GSM*geo.txt
do
	echo $sample
	awk '{printf "\t%d\t%s\n",$2,$1}' $sample > $sample.reordered.txt
done

mv GSM1551173_A671_geo.txt.reordered.txt Ctr14_GSM1551173_A671.geo.txt
mv GSM1551174_A672_geo.txt.reordered.txt Ctr16_GSM1551174_A672.geo.txt
mv GSM1551175_A673_geo.txt.reordered.txt AD6_GSM1551175_A673.geo.txt
mv GSM1551176_A675_geo.txt.reordered.txt AD5_GSM1551176_A675.geo.txt
mv GSM1551177_A676_geo.txt.reordered.txt AD3_GSM1551177_A676.geo.txt
mv GSM1551178_A677_geo.txt.reordered.txt AD1_GSM1551178_A677.geo.txt
mv GSM1551179_A678_geo.txt.reordered.txt AD4_GSM1551179_A678.geo.txt
mv GSM1551180_A679_geo.txt.reordered.txt TPD8_GSM1551180_A679.geo.txt
mv GSM1551181_A681_geo.txt.reordered.txt AD2_GSM1551181_A681.geo.txt
mv GSM1551182_A682_geo.txt.reordered.txt TPD7_GSM1551182_A682.geo.txt
mv GSM1551183_A683_geo.txt.reordered.txt Ctr15_GSM1551183_A683.geo.txt
mv GSM1551184_A684_geo.txt.reordered.txt Ctr13_GSM1551184_A684.geo.txt
mv GSM1551185_A685_geo.txt.reordered.txt Ctr12_GSM1551185_A685.geo.txt
mv GSM1551186_A686_geo.txt.reordered.txt TPD9_GSM1551186_A686.geo.txt
mv GSM1551187_A687_geo.txt.reordered.txt Ctr10_GSM1551187_A687.geo.txt
mv GSM1551188_A689_geo.txt.reordered.txt Ctr11_GSM1551188_A689.geo.txt

/projects/b1069/smallRNAscripts/combineCountsNotNorm.pl geo.txt 0 > GSE63501.publishedRawCounts.minSum0.txt
