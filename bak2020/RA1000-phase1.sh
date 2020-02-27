#########################################################################################################
######################  GWAS Study for Rheumatoid Arthritis  phase1-1000 samples #########################
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