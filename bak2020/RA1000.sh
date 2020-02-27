#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis  phase2-500 samples #########################
#########################################################################################################
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/TraylorM_31596875_GCST008993/RA_XRayDamage_European.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/TraylorM_31596875_GCST008992/RA_XRayDamage_Transethnic.txt
wget ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/OkadaY_24390342_GCST002318/RA_GWASmeta_TransEthnic_v2.txt.gz

Genetics of rheumatoid arthritis contributes to biology and drug discovery.
High-density genetic mapping identifies new susceptibility loci for rheumatoid arthritis.
Genome-wide association study meta-analysis identifies seven new rheumatoid arthritis risk loci.
Genetic associations with radiological damage in rheumatoid arthritis: Meta-analysis of seven genome-wide association studies of 2,775 cases.
Genome-wide association study in Turkish and Iranian populations identify rare familial Mediterranean fever gene (MEFV) polymorphisms associated with ankylosing spondylitis.
Genetic variation at the glycosaminoglycan metabolism pathway contributes to the risk of psoriatic arthritis but not psoriasis.
#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis  phase2-500 samples #########################
#########################################################################################################
## 02/04/2020
cd ~/hpc/rheumatology/RA/RA1000
plink --file result_extract_forward_all --make-bed --out RA1000
plink --bfile RA1000 --mind 0.02 --geno 0.1 --maf 0.01 --hwe 0.00001 --threads 31 --make-bed --out RA2020-B1
plink2 --bfile RA2020-B1 --king-cutoff 0.125
plink2 --bfile RA2020-B1 --remove plink2.king.cutoff.out.id --make-bed -out RA2020-B2
plink --bfile RA1000-B1 --maf 0.01 --make-bed --indep 50 5 2
plink --bfile RA2020-B1 --extract plink.prune.in --genome
plink --bfile RA2020-B2 --check-sex
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
head sexcheck.remove
plink --bfile RA2020-B2 --remove sexcheck.remove --make-bed --out RA2020-B3
plink --bfile RA2020-B3 --alleleACGT --snps-only just-acgt --exclude plink.dupvar --make-bed --out RA2020-B4

wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA1000/phen.pl -O phen.pl
perl phen.pl RA2020-B4.fam > RA2020-B4.fam.new
mv RA2020-B4.fam.new RA2020-B4.fam
plink --bfile RA2020-B4 --freq
plink --bfile RA2020-B4 --exclude missing.imblance.remove
plink --bfile RA2020-B4 --test-missing midp 

#### plink2vcf and Michigan imputation
mkdir michigan
cd michigan
mkdir temp
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
for i in {1..23} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile ../RA2020-B5 --chr $i --recode vcf-iid --out RA2020-B5.chr$i >> $i.job
echo bcftools view RA2020-B5.chr$i.vcf -Oz -o RA2020-B5.chr$i.vcf.gz >>$i.job
echo tabix -p vcf RA2020-B5.chr$i.vcf.gz >>$i.job
echo java -jar ./conform-gt.24May16.cee.jar gt=RA2020-B5.chr$i.vcf.gz match=POS chrom=$i ref=~/hpc/db/hg19/beagle/EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz  out=RA2020-B5.chr$i.beagle >>$i.job
echo tabix -p vcf RA2020-B5.chr$i.beagle.vcf.gz >>$i.job
qsub $i.job
done

wget https://imputationserver.sph.umich.edu/share/results/35d0da07e60e925d494f748007fb2df9/chr_1.zip 
wget https://imputationserver.sph.umich.edu/share/results/91687cd691464a9f9c80c274e57013bc/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/ca9f647452c63b5fbec5af3b2d97376d/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/881338f47b08bd96f497c8378a5bef49/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/cfa35f73333b8fd4fcfde0d0a4c2dc4b/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/223c2513ba112e1718a36280af174ae6/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/e56a9f9bd3206413af9ef27d67049325/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/54acacc4e7ddc9dd7bbefaf085f4000f/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/e67c732f70fb371f6c4cb50bedfbc744/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/bbafa82f64309b929a1385e3ae7fd026/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/9a7a5784401c26599ee7d738071fec05/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/c9dd25085f3660641a16bae9e35f32c7/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/3b43d390a2da1deaad1787ec731f48dc/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/3d60aa8084ccc306e1a2b24f338c65ac/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/5204a580d1d0a52185203ebbe8a72594/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/e4a618bcb792205a4f0172657477a113/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/bd953690afbe03f9125baeb2d9d04903/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/a98a13eaec5f26dea852d70077c05adc/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/e60c338504a56dbd9c4747beacfa92d9/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/f3814aab9d59c2602e305f54fe56861a/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/fd2b94fc77120db2ea9b79d6df5ba9b8/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/58794960f9ee17b237f00b7864f31886/chr_9.zip


cd ~/hpc/rheumatology/RA/RA1000/michigan
mkdir temp
for i in {1..22} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo unzip -P Wjq1iBSrUmx1Fp chr_$i.zip  >> $i.job
qsub $i.job
done

###### replace michigan imputation id to rs ID
mkdir temp
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
for i in {1..23} 
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=12 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo tabix -p vcf chr$i.dose.vcf.gz >>$i.job
echo java -jar ./conform-gt.24May16.cee.jar gt=chr$i.dose.vcf.gz match=POS chrom=$i ref=~/hpc/db/hg19/beagle/EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz out=chr$i.beagle >>$i.job
echo tabix -p vcf chr$i.beagle.vcf.gz >>$i.job
qsub $i.job
done

##### convert michigan imputated vcf to plink directly
mkdir temp
mkdir plink
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # gunzip chr$i.info.gz >>$i.job
echo # bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' -Oz chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz >>$i.job
echo plink --vcf chr$i.dose.vcf.gz --make-bed --out ./plink/chr$i >>$i.job
echo # plink --bfile ./plink/chr$i --list-duplicate-vars --out ./plink/chr$i >>$i.job
echo # plink --bfile chr$i --exclude chr$i.dupvar --make-bed --out ./plink/chr$i >> $i.job
qsub $i.job
done


##### convert michigan imputated vcf (pos to rs id ) and then convert to plink format
mkdir temp
mkdir beagleplink
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo # gunzip chr$i.info.gz >>$i.job
echo # bcftools view -i \'R2\>0.6\|TYPED=1\|TYPED_ONLY=1\' -Oz chr$i.dose.vcf.gz -Oz -o chr$i.dose.filter.vcf.gz >>$i.job
echo plink --vcf chr$i.beagle.vcf.gz --make-bed --out ./beagleplink/chr$i >>$i.job
echo # plink --bfile ./plink/chr$i --list-duplicate-vars --out ./plink/chr$i >>$i.job
echo # plink --bfile chr$i --exclude chr$i.dupvar --make-bed --out ./plink/chr$i >> $i.job
qsub $i.job
done

#### update fam with phenotype
wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA1000/phen.pl -O phen.pl
for i in {1..22}
do
perl phen.pl chr$i.fam > chr$i.fam.new
mv chr$i.fam.new chr$i.fam
done

#### association study with plink assoc
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out chr$i.counts >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out chr$i.freq >> $i.job
qsub $i.job
done

#### association study with plink
mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile chr$i --indep-pairwise 50 5 0.2 --out pca.chr$i >> $i.job
echo plink --bfile pca.chr$i --allow-no-sex --extract plink.prune.in --make-bed --out prune.pca.chr$i >> $i.job
echo plink --bfile pca.chr$i --pca --maf 0.05 --memory 40000 --cluster --mds-plot 4 --out pca.chr$i >> $i.job
qsub $i.job
done


touch bmerge.txt
for i in {2..22}
do
echo chr$i >> bmerge.txt
done
plink --bfile chr1 --merge-list bmerge.txt --make-bed --out RA1000

ls chr*.beagle.vcf.gz > bcfmerge.txt
bcftools concat -f bcfmerge.txt --threads 22 -Oz -o RA1000.beagle.vcf.gz
plink --vcf RA1000.beagle.vcf.gz --make-bed --out RA1000.beagle


wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA1000/phen.pl -O phen.pl
perl phen.pl RA1000.beagle.fam > RA1000.beagle.fam.new
mv RA1000.beagle.fam.new RA1000.beagle.fam
touch -m RA1000.beagle.bim
touch -m RA1000.beagle.fam
touch -m RA1000.beagle.bed

plink --bfile RA1000.beagle --indep-pairwise 50 5 0.2
plink --bfile RA1000.beagle --allow-no-sex  --extract plink.prune.in --genome --make-bed --out temp
plink --bfile temp --pca --maf 0.05 --memory 40000 --cluster --mds-plot 4

plink --bfile RA1000.beagle --allow-no-sex --logistic --threads 31 --covar plink.eigenvec --covar-number 1-4 --adjust --out RA1000.beagle

##########################################################################################
##########################################################################################
plink --bfile RA500 --allow-no-sex --logistic --threads 31 --covar RA500.eigenvec --covar-number 1-4 --adjust --out RA500
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out RA500.counts
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out RA500.freq 
cd ~/hpc/db/hg19/beagle/EUR/
wget https://raw.githubusercontent.com/zhanxw/checkVCF/master/checkVCF.py
samtools faidx hg19.fa
cd /gpfs/home/guosa/hpc/project/pmrp/cytokine/hg18
python checkVCF.py -r hs37d5.fa -o SchrodiTH17_660W.chr9  SchrodiTH17_660W.chr9.vcf.gz
plink --vcf Schrodi_IL23_IL17_combined_RECAL_SNP_INDEL_variants.VA.chr19.Minimac4.vcf.gz --double-id --freq --make-bed --out Schrodi_IL23_IL17.chr19
wget https://faculty.washington.edu/browning/beagle/beagle.16May19.351.jar -O beagle.16May19.351.jar
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar -O conform-gt.24May16.cee.jar
java -jar ./conform-gt.24May16.cee.jar gt=SchrodiTH17_660W.chr22.vcf.gz chrom=22 ref=~/hpc/db/hg19/beagle/EUR/chr22.1kg.phase3.v5a.EUR.vcf.gz  out=SchrodiTH17_660W.chr22.beagle.vcf.gz
#############################################
cd ~/hpc/db/hg19/beagle
for i in {1..22} X Y
do
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr$i.1kg.phase3.v5a.vcf.gz
done
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/20140625_related_individuals.txt
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_male_samples_v3.20130502.ALL.panel
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel
mkdir EUR
mkdir EAS

grep EUR integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}'> EUR.List.txt
grep EAS integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1}' > EAS.List.txt

mkdir temp
for i in {1..22} X Y
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo tabix -p vcf chr$i.1kg.phase3.v5a.vcf.gz >> $i.job
# echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EUR.List.txt -Oz -o ./EUR/chr$i.1kg.phase3.v5a.EUR.vcf.gz >>$i.job
echo bcftools view chr$i.1kg.phase3.v5a.vcf.gz -S EAS.List.txt -Oz -o ./EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz >>$i.job
echo bcftools norm -d all -m-both ./EAS/chr$i.1kg.phase3.v5a.dedup.EAS.vcf.gz -Oz -o ./EAS/chr$i.1kg.phase3.v5a.dedup.norm.EAS.vcf.gz  >>$i.job
qsub $i.job
done


### revise the phentoypes with R script
data1<-read.table("/gpfs/home/guosa/hpc/rheumatology/RA/RA500/RA2020-B9.fam")
file=list.files(pattern="*.fam")
for(i in file){
data2<-read.table(i)
data2[,6]<-data1[match(data2[,1],data1[,2]),6]
out=paste(i,"new",sep=".")
write.table(data2,file=out,sep=" ",quote=F,col.names=F,row.names=F)
cmd= paste("mv",out,i,sep=" ")
print(cmd)
system(cmd)
}

mkdir temp
for i in {1..22}
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=1 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out chr$i.counts >> $i.job
echo plink --bfile chr$i --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out chr$i.freq >> $i.job
qsub $i.job
done


wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/high-LD-regions.txt

plink --bfile chr1 --allow-no-sex --maf 0.05 --merge-list merge.txt --make-bed --out RA500

plink --bfile RA500 --exclude high-LD-regions.txt --range --indep-pairwise 50 5 0.2
plink --bfile RA500 --allow-no-sex  --extract plink.prune.in --genome --make-bed --out temp
plink --bfile temp --pca --maf 0.05 --memory 40000 --cluster --mds-plot 4 --out RA500
plink --bfile RA500 --allow-no-sex --logistic --threads 31 --covar RA500.eigenvec --covar-number 1-4 --adjust --out RA500
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out RA500.counts
plink --bfile RA500 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out RA500.freq 

for i in {2..22}
do
echo chr$i >> merge.txt
done

cd /gpfs/home/guosa/hpc/rheumatology/RA/RA500/michigan
gunzip chr22.info.gz
plink --vcf chr22.dose.vcf.gz --make-bed --out chr22
plink --bfile chr22 --list-duplicate-vars --out chr22
plink --bfile chr22 --exclude plink.dupvar --make-bed --out ../plink/chr22

plink --bfile RA2020 --mind 0.05 --make-bed --out RA2020-B1
plink --bfile RA2020-B1 --geno 0.95 --make-bed --out RA2020-B2
plink --bfile RA2020-B2 --maf 0.01 --make-bed --out RA2020-B3
plink --bfile RA2020-B3 --hwe 0.00001 --make-bed --out RA2020-B4
plink2 --bfile RA2020-B4 --king-cutoff 0.125
plink2 --bfile RA2020-B4 --remove plink2.king.cutoff.out.id --make-bed -out RA2020-B5
plink --bfile RA2020-B5 --check-sex
plink --bfile RA2020-B5 --impute-sex --make-bed --out RA2020-B6
plink --bfile RA2020-B6 --check-sex
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
plink --bfile RA2020-B6 --remove sexcheck.remove --make-bed --out RA2020-B7
plink --bfile RA2020-B7 --test-missing midp 
awk '$5<0.000001{print}' plink.missing | awk '{print $2}' > missing.imblance.remove
plink --bfile RA2020-B7 --exclude missing.imblance.remove --make-bed --out RA2020-B8
plink --bfile RA2020-B8 --assoc mperm=5000 --adjust gc --threads 31
plink --bfile RA2020-B8 --pca --threads 31
plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-20 --adjust
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200121.zip
##########################################################################################
##########################################################################################
plink --bfile RA2020-B4 --test-missing midp 
awk '$5<0.000001{print}' plink.missing | awk '{print $2}' > missing.imblance.remove
head missing.imblance.remove
plink --bfile RA2020-B4 --exclude missing.imblance.remove --make-bed --out RA2020-B8
plink --bfile RA2020-B8 --pca --threads 31

plink --bfile RA2020-B8 --list-duplicate-vars ids-only suppress-first
plink --bfile RA2020-B8 --maf 0.01 --make-bed --indep 50 5 2
plink --bfile RA2020-B8 --extract plink.prune.in --genome
sort -k10,10n plink.genome > plink.genome.sort

data<-read.table("plink.genome",head=T)
data<-data[order(data$PI_HAT,decreasing=T),]
head(data)
out<-data[data$PI_HAT>0.2,]
relatives<-names(tail(sort(table(c(as.character(out$IID1),as.character(out$IID2)))),10))
x1<-unique(data[which(data$IID1 %in% relatives),c(1,2)])
x2<-unique(data[which(data$IID2 %in% relatives),c(3,4)])
colnames(x1)<-c("FID","IID")
colnames(x2)<-c("FID","IID")
x<-unique(rbind(x1,x2))
write.table(x,file="plink.remove.txt",sep=" ",quote=F,col.names=F,row.names=F)
write.csv(out,file="RA1000.IBD.csv",quote=F)

tail(sort(table(c(as.character(out$IID1),as.character(out$IID2)))),10)

plink --bfile RA1000 --remove outliers.txt --make-bed --out mydata2


plink --bfile RA2020-B8 --pca --maf 0.05 --memory 40000 --cluster --mds-plot 4

wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA1000/phen.pl
perl phen.pl RA2020-B8.fam > RA2020-B8.fam.new
mv RA2020-B8.fam.new RA2020-B8.fam

plink --bfile RA2020-B8 --logistic --covar plink.eigenvec --covar-number 1-5 --adjust
plink --bfile RA2020-B8 --assoc --adjust gc --threads 31  --ci 0.95 --out RA500
plink --bfile RA2020-B8 --assoc mperm=1000000 --adjust gc --threads 31
grep "ADD\|NMISS" plink.assoc.logistic > plink.assoc.logistic.add
wget https://raw.githubusercontent.com/Shicheng-Guo/ASA/master/manhattan.plot.R -O manhattan.plot.R
Rscript manhattan.plot.R plink.assoc.logistic.add

## local
head -n 50 plink.assoc.logistic.adjusted


assoc1="/gpfs/home/guosa/hpc/rheumatology/RA/RA500/RA500.assoc"
assoc2="/gpfs/home/guosa/hpc/rheumatology/RA/he2020/RA1000.assoc"
plink --meta-analysis $assoc1 $assoc2

cd ../he2020
plink --bfile RA2020-B8 --assoc --adjust gc --threads 31  --ci 0.95  --counts --out RA1000.counts
plink --bfile RA2020-B8 --assoc --adjust gc --threads 31  --ci 0.95  --out RA1000

cd ../RA500
plink --bfile RA2020-B8 --assoc --adjust gc --threads 31  --ci 0.95  --counts --out RA500.counts
plink --bfile RA2020-B8 --assoc --adjust gc --threads 31  --ci 0.95  --out RA500


/gpfs/home/guosa/hpc/tools/RSEM_tutorial/software/RSEM-1.2.25

