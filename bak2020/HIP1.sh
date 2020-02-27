
cd ~/hpc/rheumatology/RA/meta3000
input="HIP1"
mkdir $input
cd $input
grep -w "HIP1" ~/hpc/db/hg19/refGene.hg19.V2.bed.txt | awk '{print $1,$2-50000,$3+50000,$5}' OFS="\t" | sort -u | bedtools sort -i > ROI.hg19.bed

# awk '{print $3,$5,$6,$4,$2,$13}' OFS="\t" refGene.txt > ~/hpc/db/hg19/refGene.hg19.fancylocal.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/localhit.pl -O localhit.pl
wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/R/make.fancy.locus.plot.unix.R -O make.fancy.locus.plot.unix.R

####################################################################################################################
####################################################################################################################
# RA500
plink --bfile ~/hpc/rheumatology/RA/RA500/michigan/R500.beagle --extract range ROI.hg19.bed --make-bed --out ROI.RA500.dbsnp
plink --bfile ROI.RA500.dbsnp --allow-no-sex --hardy --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --assoc --counts --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --fisher --counts --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --model fisher --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA500.dominant
plink --bfile ROI.RA500.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA500.recessive
# ~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA500.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA500.haplotype --noweb
awk '$6=="NMISS"{print}' ROI.RA500.assoc.logistic > ROI.RA500.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.RA500.assoc.logistic >> ROI.RA500.assoc.logistic.add
grep -v NA ROI.RA500.assoc.logistic | sort -k12n,12 | head -n 5
grep -v NA ROI.RA500.assoc.logistic | sort -k12n,12 | tail -n 5
####################################################################################################################
####################################################################################################################
# RA1000
plink --bfile ~/hpc/rheumatology/RA/RA1000/michigan/RA1000.beagle --extract range ROI.hg19.bed --make-bed --out ROI.RA1000.dbsnp
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --hardy --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --assoc --counts --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --fisher --counts --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --model fisher --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA1000.dominant
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA1000.recessive
# ~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA1000.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA1000.haplotype --noweb
awk '$6=="NMISS"{print}' ROI.RA1000.assoc.logistic > ROI.RA1000.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.RA1000.assoc.logistic >> ROI.RA1000.assoc.logistic.add
grep -v NA ROI.RA1000.assoc.logistic | sort -k12n,12 | head -n 5
grep -v NA ROI.RA1000.assoc.logistic | sort -k12n,12 | tail -n 5
####################################################################################################################
####################################################################################################################
# plink --bfile ~/hpc/rheumatology/RA/RA1000/michigan/RA1000.beagle --bmerge ~/hpc/rheumatology/RA/RA500/michigan/R500.beagle --make-bed --out ~/hpc/rheumatology/RA/meta300/RA3000.beagle 
# plink --bfile ~/hpc/rheumatology/RA/RA1000/michigan/RA1000.beagle --exclude RA3000.beagle-merge.missnp --make-bed --out RA1000.beagle.demiss
# plink --bfile ~/hpc/rheumatology/RA/RA500/michigan/R500.beagle --exclude RA3000.beagle-merge.missnp --make-bed --out RA500.beagle.demiss
# plink --bfile RA1000.beagle.demiss --bmerge RA500.beagle.demiss --allow-no-sex --make-bed --out RA3000.beagle 
####################################################################################################################
####################################################################################################################
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --extract range ROI.hg19.bed --make-bed --out ROI.RA3000.dbsnp
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --hardy --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --assoc --counts --adjust --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --fisher --counts --adjust --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --model fisher --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA3000.dominant
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA3000.recessive
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA3000.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA3000.haplotype --noweb
awk '$6=="NMISS"{print}' ROI.RA3000.assoc.logistic > ROI.RA3000.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.RA3000.assoc.logistic >> ROI.RA3000.assoc.logistic.add
grep -v NA ROI.RA3000.assoc.logistic | sort -k12n,12 | head -n 5
grep -v NA ROI.RA3000.assoc.logistic | sort -k12n,12 | tail -n 5

rs="rs140115892"
chr=4
grep $rs ROI.RA3000.assoc.logistic
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --maf 0.01 --allow-no-sex --snp $rs --recode --out $rs
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --maf 0.01 --allow-no-sex --threads 31 --r2 --ld-snp $rs --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0 --out $rs
perl localhit.pl $rs > $rs.local
head  $rs*
Rscript make.fancy.locus.plot.unix.R $rs $rs $chr $rs.local 6 0.05
