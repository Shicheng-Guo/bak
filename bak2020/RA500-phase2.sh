#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis  phase2-500 samples #########################
#########################################################################################################
## 02/04/2020
cd ~/hpc/rheumatology/RA/RA500
scp -o 'ProxyCommand ssh nu_guos@submit-1.chtc.wisc.edu nc %h %p' root@101.133.145.142:/data/ASA-CHIA_20191211/second_time/data_result/result* ./
plink --bfile result_extract_forward --make-bed --out RA500
plink --bfile RA500 --mind 0.02 --geno 0.1 --maf 0.01 --hwe 0.00001 --threads 31 --make-bed --out RA500-B1
plink2 --bfile RA500-B1 --king-cutoff 0.125
plink --bfile RA500-B1 --genome --min 0.25
plink2 --bfile RA500-B1 --remove plink2.king.cutoff.out.id --make-bed -out RA500-B2
plink --bfile RA500-B1 --maf 0.01 --make-bed --indep 50 5 2
plink --bfile RA500-B1 --extract plink.prune.in --genome
plink --bfile RA500-B2 --check-sex
grep PROBLEM plink.sexcheck 
grep PROBLEM plink.sexcheck | awk '{print $1,$2}' > sexcheck.remove
plink --bfile RA500-B2 --remove sexcheck.remove --make-bed --out RA500-B3
plink --bfile RA500-B3 --list-duplicate-vars ids-only suppress-first
plink --bfile RA500-B3 --alleleACGT --snps-only just-acgt --exclude plink.dupvar --make-bed --out RA500-B4

wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA500/phen.pl -O phen.pl
perl phen.pl RA500-B4.fam > RA500-B4.fam.new
mv RA500-B4.fam.new RA500-B4.fam
plink --bfile RA500-B4 --freq
plink --bfile RA500-B4 --test-missing midp 
plink --bfile RA500-B4 --exclude missing.imblance.remove

plink --bfile RA500-B4 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out RA500-B4.counts
plink --bfile RA500-B4 --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out RA500-B4.freq
plink --bfile R500.beagle --allow-no-sex --logistic --threads 31 --covar plink.eigenvec --covar-number 1-4 --adjust --out R500.beagle

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
echo plink --bfile ../RA500-B4 --chr $i --recode vcf-iid --out RA500-B4.chr$i >> $i.job
echo bcftools view RA500-B4.chr$i.vcf -Oz -o RA500-B4.chr$i.vcf.gz >>$i.job
echo tabix -p vcf RA500-B4.chr$i.vcf.gz >>$i.job
echo java -jar ./conform-gt.24May16.cee.jar gt=RA500-B4.chr$i.vcf.gz match=POS chrom=$i ref=~/hpc/db/hg19/beagle/EAS/chr$i.1kg.phase3.v5a.EAS.vcf.gz  out=RA500-B4.chr$i.beagle >>$i.job
echo tabix -p vcf RA500-B4.chr$i.beagle.vcf.gz >>$i.job
qsub $i.job
done

wget https://imputationserver.sph.umich.edu/share/results/14b1257e930de8629ee101a149050b29/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/1b77872f088f4dcfa12ba793fcd7c38e/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/d9808557b31f531a98892c3729b71a65/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/f64d326014a15308048cde335fc3804a/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/1b0b93646ca044ee86f5c75af7954fe7/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/2254a1c8bca44c0fd5e7e4addf2cc34e/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/75a5a53715bf760d6750fae28b992cfc/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/e5e546dcd69b7ec3bfaae08d06e1e8b9/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/a80c04ef078409539a6cc7e1d48f52af/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/179a27961700c8bce5692d557d2c4b40/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/9f4b56292b7b7a42707df35b2a76d8fc/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/6b48804237ff92c12c8c3b2c085a1205/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/eb0c063eeecb1fbb0b0427275adadbea/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/a6cf6679fe5d5aacf33eafc18d3d015c/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/b11b7b175d53e7ca37d1bca24cbcfb41/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/2f6b559c284adacce44467a273964d8d/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/1281d92bdc10700719999bb3460563e2/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/c02f0976267ccf6b772975c94c650aad/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/8befd8d955c122a3c5ec5c56edfaeb41/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/2a2b56fdc5456e3e2d48cd949c1c23cd/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/e5be2b1d92c5927dade556a84f803bab/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/71c5777e6e2a862775d3ddbc0a167aa9/chr_9.zip

cd ~/hpc/rheumatology/RA/RA500/michigan
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
echo unzip -P pom36AZlvIUYC8 chr_$i.zip  >> $i.job
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

## merge vcf with chrosome
ls chr*.beagle.vcf.gz > bcfmerge.txt
bcftools concat -f bcfmerge.txt --threads 32 -Oz -o R500.beagle.vcf.gz
plink --vcf R500.beagle.vcf.gz --make-bed --out R500.beagle

wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/RA500/phen.pl -O phen.pl
perl phen.pl R500.beagle.fam > R500.beagle.fam.new
head R500.beagle.fam.new

mv R500.beagle.fam.new R500.beagle.fam
touch -m R500.beagle.bim
touch -m R500.beagle.fam
touch -m R500.beagle.bed

plink --bfile R500.beagle --indep-pairwise 50 5 0.2
plink --bfile R500.beagle --allow-no-sex  --extract plink.prune.in --genome --make-bed --out temp
plink --bfile temp --pca --threads 31 --maf 0.05 --memory 40000 --cluster --mds-plot 4

plink --bfile R500.beagle --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --counts --out R500.beagle.counts 
plink --bfile R500.beagle --assoc --adjust gc --threads 31 --allow-no-sex  --ci 0.95  --out R500.beagle.freq 
plink --bfile R500.beagle --allow-no-sex --logistic --threads 31 --covar plink.eigenvec --ci 0.95  --covar-number 1-4 --adjust --out RA500.beagle

############################################################################################################################################################
############################################################################################################################################################
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

