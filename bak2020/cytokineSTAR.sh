## build STAR Genome Index 
wget http://ftp.ensemblorg.ebi.ac.uk/pub/grch37/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.transcripts.fa.gz
gunzip gencode.v32lift37.transcripts.fa.gz
gunzip Homo_sapiens.GRCh37.87.gtf.gz
mv gencode.v32lift37.annotation.gtf ~/hpc/db/hg19/STAR/
##
STAR --runMode genomeGenerate --genomeDir ~/hpc/db/hg19/STAR/ --genomeFastaFiles ~/hpc/db/hg19/hg19.fa --sjdbGTFfile ~/hpc/db/hg19/gencode.v29lift37.annotation.gtf --runThreadN 30 --sjdbOverhang 89
echo chrLength.txt
echo chrNameLength.txt
echo chrName.txt
echo chrStart.txt
## Mapping RNA-seq with STAR pipeline
cd /gpfs/home/guosa/hpc/project/pmrp/cytokine/rnaseq/bam
mkdir  AC90P8ANXX
mkdir AC90KJANXX
mkdir AC907MANXX
cd /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/AC90P8ANXX
cd /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/AC907MANXX
cd /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/AC90KJANXX

mkdir temp
OPTS_P1="--outSAMstrandField intronMotif --twopassMode Basic --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000 --genomeLoad NoSharedMemory --seedSearchStartLmax 8 --outFilterMismatchNoverLmax 0.5"
OPTS_P2="--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000" 
DBDIR="~/hpc/db/hg19/STAR"
ZCAT="--readFilesCommand zcat"
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
#echo bowtie2 -p 6 -x /gpfs/home/guosa/hpc/db/hg19/bowtie2/hg19 -1 $i\_1.fastq.gz -2 $i\_2.fastq.gz -S $i.sam >> $i.job
#echo samtools view -bS $i.sam \> $i.bam >> $i.job
#echo samtools sort $i.bam -o $i.sorted.bam >> $i.job
#echo samtools mpileup -uf ~/hpc/db/hg19/hg19.db $i.sorted.bam \| bcftools view -Ov - \> $i.bcf >> $i.job
#echo samtools depth $i.sorted.bam \> $i.wig >> $i.job
# echo STAR $OPTS_P1 $OPTS_P2 --runThreadN 24 --outBAMsortingThreadN 6 $ZCAT --genomeDir $DBDIR  --sjdbGTFfile ~/hpc/db/hg19/STAR/gencode.v32lift37.annotation.gtf --outFileNamePrefix ~/hpc/project/pmrp/cytokine/rnaseq/bam/AC90KJANXX/$i --readFilesIn $i\_R1_001.fastq.gz $i\_R2_001.fastq.gz >> $i.job
#echo mkdir ~/hpc/project/pmrp/cytokine/rnaseq/bam/cufflinks/$i >> $i.job
#echo cufflinks -p 12 --GTF ~/hpc/db/hg19/STAR/gencode.v32lift37.annotation.gtf ~/hpc/project/pmrp/cytokine/rnaseq/bam/AC90P8ANXX/$i\Aligned.sortedByCoord.out.bam -o ~/hpc/project/pmrp/cytokine/rnaseq/bam/cufflinks/$i >> $i.job
# echo trim_galore --paired --phred33 --clip_R1 2 --clip_R2 2 --fastqc --illumina $i\_R1_001.fastq.gz $i\_R2_001.fastq.gz --output_dir ./fastq_trim >> $i.job
# echo fastq_screen --subset 1000000 --force --threads 12 --conf  ~/hpc/tools/fastq_screen_v0.14.0/FastQ_Screen_Genomes/fastq_screen.conf $i\_R1_001.fastq.gz $i\_R2_001.fastq.gz --outdir ./fastqscreen/ >> $i.job
qsub  $i.job
done

UW040LPS_ATTCCT_L002Aligned.sortedByCoord.out.bam

htseq-count -m intersection-nonempty -t exon -i gene_id -f bam starGC025680Aligned.toTranscriptome.out.bam ~/data/genomes/pig/susScr3/star/susScr3.04012016.gtf -o Test

cd ~/hpc/project/pmrp/cytokine/rnaseq/bam/AC90P8ANXX/
cd ~/hpc/project/pmrp/cytokine/rnaseq/bam/AC907MANXX/
cd ~/hpc/project/pmrp/cytokine/rnaseq/bam/AC90KJANXX/

cd /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/AC90P8ANXX

salmon quant -t transcripts.fa -l A -a UW040LPS_ATTCCT_L002Aligned.sortedByCoord.out.bam -o salmon_quant

infer_experiment.py -r ../knowngene.hg19.bed12 -i UW040LPS_ATTCCT_L002Aligned.sortedByCoord.out.bam


cd AC90P8ANXX
rm ../mappingratio.csv
grep 'Uniquely mapped reads %' *.final.out | awk -F" " '{print $1,$2,$7}' OFS="," > ../mappingratio.csv
cd ../AC907MANXX
grep 'Uniquely mapped reads %' *.final.out | awk -F" " '{print $1,$2,$7}' OFS="," >> ../mappingratio.csv
cd ../AC90KJANXX
grep 'Uniquely mapped reads %' *.final.out | awk -F" " '{print $1,$2,$7}' OFS="," >> ../mappingratio.csv


cd ~/hpc/project/pmrp/cytokine/rnaseq/

cp /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/*/*fastq.gz ./

cd /gpfs/home/guosa/hpc/project/pmrp/cytokine/rnaseq

cd /mnt/bigdata/Genetic/Projects/Schrodi_IL23_IL17_variants/RNAseq_macrophages/Data/

/gpfs/home/guosa/hpc/project/pmrp/cytokine/rnaseq

484

cd AC90P8ANXX
ls *.fastq.gz > ../filename.txt
cd ../AC90KJANXX
ls *.fastq.gz >> ../filename.txt
cd ../AC907MANXX
ls *.fastq.gz >> ../filename.txt

sort -u filename.txt | wc -l 




ERROR: deeptools 3.2.0 requires deeptoolsintervals>=0.1.7, which is not installed.
ERROR: deeptools 3.2.0 requires matplotlib>=2.1.2, which is not installed.
ERROR: deeptools 3.2.0 requires numpydoc>=0.5, which is not installed.
ERROR: deeptools 3.2.0 requires plotly>=2.0.0, which is not installed.
ERROR: deeptools 3.2.0 requires py2bit>=0.2.0, which is not installed.
ERROR: deeptools 3.2.0 requires pyBigWig>=0.2.1, which is not installed.
ERROR: deeptools 3.2.0 requires scipy>=0.17.0, which is not installed.

pip install deeptoolsintervals
pip install matplotlib
pip install numpydoc
pip install plotly
pip install py2bit
pip install pyBigWig
pip install scipy

