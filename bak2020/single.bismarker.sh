#!/bin/csh
#PBS -N SRR10775514
#PBS -q shortq
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
#PBS -o SRR10775514.log
#PBS -e SRR10775514.err
#PBS -V
#PBS -M shihcheng.guo@gmail.com
#PBS -m abe
#PBS -A k4zhang-group
cd /gpfs/home/guosa/PRJNA597894/sra
fastq-dump --skip-technical --gzip SRR10775514
trim_galore --phred33 -e 0.2 --fastqc SRR10775514.fastq.gz --output_dir ../fastq_trim
bismark --bowtie2 --phred33-quals --fastq -N 1 --multicore 2 ~/hpc/db/hg19/bismark ../fastq_trim/SRR10775514_trimmed.fq.gz -o ../bam
filter_non_conversion --single ../bam/SRR10775514_trimmed_bismark_bt2.bam
deduplicate_bismark --bam ../bam/SRR10775514_trimmed_bismark_bt2.nonCG_filtered.bam
bismark_methylation_extractor --single-end --merge_non_CpG --bedGraph --cutoff 5 --ignore 1 --buffer_size 4G --comprehensive --output ../methyfreq  ../bam/SRR10775514_trimmed_bismark_bt2.nonCG_filtered.deduplicated.bam
samtools sort ../bam/SRR10775514_trimmed_bismark_bt2.nonCG_filtered.deduplicated.bam -o ../sortbam/SRR10775514.bismark_bt2_se.sort.bam
samtools index ../sortbam/SRR10775514.bismark_bt2_se.sort.bam
perl ~/bin/samInfoPrep4Bam2Hapinfo.pl ~/oasis/db/hg19/hg19.cut10k.bed > saminfo.txt
perl ~/bin/bam2hapInfo2PBS.pl saminfo.txt submit bismark ~/hpc/db/hg19/hg19.chrom.sizes ~/hpc/db/hg19/HsGenome19.CpG.positions.txt
