#########################################################################################################
##### RNA-seq data to reveal novel response mechanism to bacterial within host wound tissues ###########
#########################################################################################################
## 02/12/2020
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("d3heatmap")
library("dendextend")
library("ComplexHeatmap")
library("circlize")

for i in $(ls *.fastq | rev | cut -c 19- | rev | uniq)
do
mkdir $i
mv $i* ./$i
done

for i in $(cat xx)
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo cd $(pwd) >> $i.job
echo sh \.\/samsa2.sh $i >>$i.job
qsub  $i.job
done



for i in `ls *.fastq`
do
echo $i
echo \#PBS -N $i > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo cd $(pwd) >> $i.job
echo metaphlan2.py $i --input_type fastq --nproc 24 \> $i\_profile.txt >>$i.job
# qsub  $i.job
done

#########################################################################################################
##### RNA-seq data to reveal novel response mechanism to bacterial within host wound tissues ###########
#########################################################################################################
## 02/04/2020
cd ~/hpc/project/RnaseqBacterial/extdata/rnaseq
wget http://cs.wellesley.edu/~btjaden/Rockhopper/download/current/Rockhopper.jar
cd ~/hpc/project/RnaseqBacterial/extdata/rnaseq
genome_DIR1=~/hpc/project/RnaseqBacterial/extdata/rnaseq/Rockhopper_Results/genomes/Staphylococcus_aureus_subsp__aureus_USA300_FPR3757
mkdir temp
mkdir diamond
for i in $(ls *.fastq.gz | rev | cut -c 17- | rev | uniq)
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo java -Xmx1200m -cp Rockhopper.jar Rockhopper -g $genome_DIR1 $i\_R1_001.fastq.gz%$i\_R2_001.fastq.gz -o ./Rockhopper_Results/$i >> $i.job
echo pear -f $i\_R1_001.fastq.gz -r $i\_R2_001.fastq.gz -o $i >> $i.job
echo zcat $i\_R1_001.fastq.gz \> $i.fastq >> $i.job
echo zcat $i\_R2_001.fastq.gz \>\> $i.fastq >> $i.job
echo diamond blastx -d ~/hpc/db/nr -q $i.fastq -o ./diamond/$i >> $i.job
qsub  $i.job
done

cp *L001.fastq /gpfs/home/guosa/hpc/tools/samsa2/input_files 
cp *_001.fastq.gz /gpfs/home/guosa/hpc/tools/samsa2/input_files 

python ~/hpc/tools/samsa2/python_scripts/DIAMOND_analysis_counter.py -I 42_S3_L001 -D ~/hpc/db/nr -O 

PATH=~/hpc/tools/FastQC:$PATH
PATH=~/hpc/tools/RSEM_tutorial/software/RSEM-1.2.25:$PATH
PATH=~/hpc/tools/sortmerna-2.1:$PATH
PATH=~/hpc/tools/pear-0.9.11-linux-x86_64/bin:$PATH
PATH=~/hpc/tools/TrimGalore-0.4.5:$PATH
#########################################################################################################
##### RNA-seq data to reveal novel response mechanism to bacterial within host wound tissues ###########
#########################################################################################################
## 02/04/2020

mkdir ~/hpc/tools/diamond
cd ~/hpc/tools/diamond

wget http://github.com/bbuchfink/diamond/releases/download/v0.9.30/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr

git clone https://github.com/transcript/samsa2.git
cd setup_and_test
sh ./package_installation.bash

wget https://github.com/Shicheng-Guo/RnaseqBacterial/blob/master/bin/pear-0.9.11-linux-x86_64.tar.gz
tar xzvf pear-0.9.11-linux-x86_64.tar.gz

wget https://github.com/biocore/sortmerna/archive/2.1.tar.gz

diamond makedb --in nr.faa -d nr
diamond blastx -d nr -q reads.fna -o matches.m8

#########################################################################################################
##### SAMSA2:  ###########
#########################################################################################################
## 02/04/2020
cd ~/hpc/tools/samsa2
wget https://raw.githubusercontent.com/Shicheng-Guo/RnaseqBacterial/master/bin/master_script_preserving_unmerged.sh

cd ~/hpc/tools/samsa2/bash_scripts
rm ~/hpc/tools/samsa2/input1/checkpoints
rm -rf ~/hpc/tools/samsa2/result1/*
~/hpc/tools/samsa2/bash_scripts/master_script_preserving_unmerged.sh
PEAR v0.9.10 was applied to mrege paired-end Rnaseq reads and merged with non-merged reads as the single-end sequencing dataset. 
perl xls2matrix.pl > strain.tab.matrix.freq.txt




vim .pl


