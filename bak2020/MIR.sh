wget https://raw.githubusercontent.com/Shicheng-Guo/miRNA-RA/master/db/hsa.gff.hg19.bed -O hsa.gff3.hg38.bed 
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O liftOver
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz -O hg18ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -O hg19ToHg38.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -O hg38ToHg19.over.chain.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/GscRbasement/master/manhattan.qqplot.R -O manhattan.plot.R
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/localhit.pl -O localhit.pl
wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/R/make.fancy.locus.plot.unix.R -O make.fancy.locus.plot.unix.R
####################################################################################################################
####################################################################################################################
cd ~/hpc/rheumatology/RA/meta3000
input="MIR"
mkdir $input
cd $input
# grep -w "MIR" ~/hpc/db/hg19/refGene.hg19.V2.bed.txt | awk '{print $1,$2,$3,$5}' OFS="\t" | sort -u | bedtools sort -i > ROI.hg19.bed
# wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/meta/ARD2019/ARD2019.hg19.txt -O ARD2019.hg19.txt 
gunzip hg18ToHg19.over.chain.gz
gunzip hg19ToHg38.over.chain.gz
gunzip hg38ToHg19.over.chain.gz
grep -v '#' hsa.gff3 | awk '{print $1,$4,$5,$9}'  OFS="\t" > hsa.gff3.hg38.bed 
liftOver hsa.gff3.hg38.bed hg38ToHg19.over.chain hsa.gff3.hg19.txt unmap
grep -v 'Un' hsa.gff3.hg19.txt > hsa.gff3.hg19.bed
perl -p -i -e 's/chr//g' hsa.gff3.hg19.bed
perl -p -i -e 's/X/23/g' hsa.gff3.hg19.bed
perl -p -i -e 's/Y/24/g' hsa.gff3.hg19.bed
wget https://raw.githubusercontent.com/Shicheng-Guo/rheumatoidarthritis/master/extdata/meta/MIR/CTRL.SNP.hg19.bed -O CTRL.SNP.hg19.bed

d1<-unique(read.table("miRNA.ESA.gnomad.common.uni.hg19.bed",sep="\t"))
d2<-unique(read.table("miRNA.ESA.GenomeAsia100K.common.hg19.bed.uni.bed",sep="\t"))
d3<-unique(d1[d1$V4 %in% d2$V4,])
write.table(d3,file="hsa.gff3.common.hg19.bed",sep="\t",quote=F,col.names=F,row.names=F)

head hsa.gff3.common.hg19.bed
wc -l hsa.gff3.common.hg19.bed
cat CTRL.SNP.hg19.bed >> hsa.gff3.common.hg19.bed
wc -l hsa.gff3.common.hg19.bed
cp hsa.gff3.common.hg19.bed hsa.gff3.hg19.bed
# wget https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/miRNA/miRNA-SNP/miRNAcommonSNP.txt -O miRNAcommonSNP.txt
# awk '{print $3}' miRNAcommonSNP.txt > miRNAcommon.hg19.txt
# awk '{print $3,$5,$6,$4,$2,$13}' OFS="\t" refGene.txt > ~/hpc/db/hg19/refGene.hg19.fancylocal.bed
####################################################################################################################
####################################################################################################################
# vcf2bed
for i in `ls gnomad.genomes.r2.1.sites.chr*.rec.hsa.gff3.sort.rmdup.biallelic.vcf.bgz`
do 
bcftools query -f '%CHROM\t%POS\t%ID\n' $i | awk '{print $1,$2-1,$2,$3}' OFS="\t"
done
bcftools view -m2 -M2 -v snps -R /home/guosa/hpc/rheumatology/RA/meta3000/MIR/hsa.gff3.hg19.bed dbSNP153.hg19.vcf.gz -Oz -o ~/hpc/rheumatology/RA/meta3000/dbSNP153.miR.hg19.vcf.gz
####################################################################################################################
####################################################################################################################
# RA500
plink --bfile ~/hpc/rheumatology/RA/RA500/michigan/R500.beagle --extract range hsa.gff3.hg19.bed --make-bed --out ROI.RA500.dbsnp
plink --bfile ROI.RA500.dbsnp --allow-no-sex --hardy --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI
plink --bfile ROI.RA500.dbsnp --allow-no-sex --assoc --counts --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --freq counts --assoc fisher --adjust --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --freq counts --model mperm=1000000 --ci 0.95 --out ROI.RA500
plink --bfile ROI.RA500.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA500.dominant
plink --bfile ROI.RA500.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA500.recessive
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA500.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA500.haplotype --noweb
awk '$6=="NMISS"{print}' ROI.RA500.assoc.logistic > ROI.RA500.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.RA500.assoc.logistic >> ROI.RA500.assoc.logistic.add
grep -v NA ROI.RA500.assoc.logistic | sort -k12n,12 | head -n 5
grep -v NA ROI.RA500.assoc.logistic | sort -k12n,12 | tail -n 5
####################################################################################################################
####################################################################################################################
# RA1000
plink --bfile ~/hpc/rheumatology/RA/RA1000/michigan/RA1000.beagle --extract range hsa.gff3.hg19.bed --make-bed --out ROI.RA1000.dbsnp
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --hardy --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --assoc --counts --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --freq counts --assoc fisher --adjust --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --freq counts --model mperm=1000000 --ci 0.95 --out ROI.RA1000
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA1000.dominant
plink --bfile ROI.RA1000.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA1000.recessive
~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA1000.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA1000.haplotype --noweb
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
# plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --extract HLA-DRB1.eqtl.txt --maf 0.01 --make-bed --out ROI.RA3000.dbsnp
####################################################################################################################
####################################################################################################################
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --extract range hsa.gff3.hg19.bed --make-bed --out ROI.RA3000.dbsnp
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --hardy --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --counts --maf 0.01 --logistic --hwe 0.01 --adjust --ci 0.95 --out ROI

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --maf 0.01 --assoc counts --hwe 0.01 --adjust --ci 0.95 --out ROI.RA3000.counts
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --maf 0.01 --assoc --adjust --ci 0.95 --out ROI.RA3000.freq


plink --bfile ROI.RA3000.dbsnp --allow-no-sex --extract MySigSNPs.txt --reference-allele MyReferenceAllele.txt --maf 0.01 --model fisher --hwe 0.01 --ci 0.95 --out ROI.RA3000

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --extract MySigSNPs.txt --reference-allele MyReferenceAllele.txt --hwe 0.01 --maf 0.01 --recode --out ROI.RA3000

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --extract MySigSNPs.txt --reference-allele MyReferenceAllele.txt --recode12 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --recodeA --out ROI.RA3000




plink2 --pfile ROI.RA3000.dbsnp --glm dominant --reference-allele MyReferenceAllele.txt --maf 0.01   --out ROI.RA3000.freq

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --maf 0.01 --reference-allele MyReferenceAllele.txt --epistasis --epi1 0.001 --out ROI.RA3000

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --freq --logistic --extract MySigSNPs.txt --model fisher --ci 0.95 --out ROI.RA3000
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --logistic --dominant --ci 0.95 --out ROI.RA3000.dominant
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --logistic --recessive --ci 0.95 --out ROI.RA3000.recessive
#~/hpc/tools/plink-1.07-x86_64/plink --bfile ROI.RA3000.dbsnp --hap-window 2,3,4,5,6 --hap-assoc --out ROI.RA3000.haplotype --noweb
awk '$6=="NMISS"{print}' ROI.RA3000.assoc.logistic > ROI.RA3000.assoc.logistic.add
awk '$5=="ADD"{print}' ROI.RA3000.assoc.logistic >> ROI.RA3000.assoc.logistic.add
grep -v NA ROI.RA3000.assoc.logistic | sort -k12n,12 | head -n 10
grep -v NA ROI.RA3000.assoc.logistic | sort -k12n,12 | tail -n 10

plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --threads 31 --epistasis --twolocus rs1414273  rs117344178  --out rs1414273.rs117344178
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --threads 31 --epistasis --twolocus rs6997249  rs2273626  --out rs6997249.rs2273626
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --threads 31 --epistasis --twolocus rs2114358  rs76468441  --out rs2114358.rs76468441
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --threads 31 --epistasis --twolocus rs4285314  rs16958290  --out  rs4285314.rs16958290
plink --bfile ROI.RA3000.dbsnp --allow-no-sex --reference-allele MyReferenceAllele.txt --threads 31 --epistasis --twolocus rs367805  rs641071  --out rs367805.rs641071

## Method 1: This method is quite slow, not recommend. 
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/rs2bed.pl -O rs2bed.pl
perl rs2bed.pl ctrl.txt
## Method 1: This method is quite slow, not recommend. 
rs="rs2620381"
chr=4
grep $rs ROI.RA3000.assoc.logistic
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --maf 0.01 --allow-no-sex --snp $rs --recode --out $rs
plink --bfile ~/hpc/rheumatology/RA/meta3000/RA3000.beagle --maf 0.01 --allow-no-sex --threads 31 --r2 --ld-snp $rs --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0 --out $rs
perl localhit.pl $rs > $rs.local
head  $rs*
Rscript make.fancy.locus.plot.unix.R $rs $rs $chr $rs.local 6 0.05

grep rs9275376 ROI.RA3000.assoc.logistic
grep rs174583 ROI.RA3000.assoc.logistic
grep rs11745587 ROI.RA3000.assoc.logistic
grep rs8005568 ROI.RA3000.assoc.logistic
grep rs6117562 ROI.RA3000.assoc.logistic
grep rs4941430 ROI.RA3000.assoc.logistic

setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/rheumatology/RA/meta3000/MIR")
file=list.files(pattern="model$")
data<-read.table(file[3],head=T,sep="")  
data$key=paste(data$SNP,data$TEST,sep=":")
for(i in 1:2){
  temp<-read.table(file[i],head=T,sep="")  
  temp$key=paste(temp$SNP,temp$TEST,sep=":")
  data<-merge(data,temp,by="key")
}
head(data)
frame<-data.frame(data$P.x,data$P.y,data$P)
colnames(data)<-gsub(".x",".RA500",colnames(data))
colnames(data)<-gsub(".y",".RA1000",colnames(data))
subset(data,P.RA500<0.05&P<0.01)

ls gnomad.genomes.r2.1.sites.chr*.rec.hsa.gff3.sort.rmdup.biallelic.vcf.bgz > xx.txt
bcftools concat -f xx.txt -Oz -o gnomad.genomes.r2.1.sites.rec.hsa.gff3.sort.rmdup.biallelic.vcf.gz

for i in `ls gnomad.genomes.r2.1.sites.chr*.rec.hsa.gff3.sort.rmdup.biallelic.vcf.bgz`
do 
bcftools query -f '%CHROM\t%POS\t%ID\n' $i | awk '{print $1,$2-1,$2,$3}' OFS="\t"
done

####################################################################################################################
####################################################################################################################
# PCA analysis
awk '{print $4}' hsa.gff3.hg19.bed  | sort -u > hsa.gff3.hg19.snp
plink --bfile /gpfs/home/guosa/hpc/db/hg19/1000Genome/plink/G1000plink --extract hsa.gff3.hg19.snp --make-bed --out G1000.miRNA
plink --bfile ROI.RA3000.dbsnp --bmerge G1000.miRNA --allow-no-sex --make-bed --out ./PCA/ROI
plink --bfile ROI --threads 31 --cluster --mds-plot 2
plink --bfile ROI --threads 31 --pca 2 'header' --out ROI

hapmap2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hapmap2.pop",head=F)
hapmap3<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hapmap3.pop",head=T)
G1000Sam <-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/1000G/1000GenomeSampleInfo.txt",head=F,as.is=T)
G1000Super<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/1000G/superpopulation.txt",head=F,sep="\t")
G1000Sam$superpop<-G1000Super[match(G1000Sam$V3,G1000Super$V1),]$V2
write.table(G1000Sam,file="1000GenomeSampleInfo.txt",quote=F,sep="\t",col.names=F,row.names=F)

eigenvec<-read.table("plink.eigenvec",head=T)
head(eigenvec)
pop<-as.character(G1000Sam[match(eigenvec[,2],G1000Sam[,2]),3])
super<-as.character(G1000Sam[match(eigenvec[,2],G1000Sam[,2]),5])
pop[is.na(pop)]<-"GHRA"
super[is.na(super)]<-"GHRA"
eigenvec$pop=pop
eigenvec$super=super
eigenvec$col=as.numeric(as.factor(super))
eigenvec$pch=as.numeric(as.factor(super))

set<-unique(data.frame(pch=eigenvec$pch,col=eigenvec$col,legend=eigenvec$super))

for(i in 1:15){
jpeg(paste("pca.super",i,".jpg",sep=""))
plot(eigenvec[,3],eigenvec[,4],pch=16,col=eigenvec$col+i,xlab="principle component 1",ylab="principle componment 2",cex.axis=1.5,cex.lab=1.5,cex=1)
legend("topright",pch=16,legend=set$legend,col=set$col+i,bty="n",cex=1)
dev.off()
}

####################################################################################################################
####################################################################################################################
# Prepare Tables and Figures with R script
BiocManager::install("SNPassoc")

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
data2<-read.table("ROI.RA3000.freq.assoc",head=T,as.is=T)
data1<-data1[order(data1$P,decreasing=F),]
data2<-data2[order(data2$P,decreasing=F),]
data2$F_A=percent(data2$F_A)
data2$F_U=percent(data2$F_U)
data2$A1=paste(data2$A1,"/",data2$A2,sep="")
data2$F_A=paste(data2$F_A,"(",data1$C_A,")",sep="")
data2$F_U=paste(data2$F_U,"(",data1$C_U,")",sep="")
data2$OR=paste(round(data2$OR,2),"(",round(data2$L95,2),"-",round(data2$U95,2),")",sep="")
data2<-data2[,-c(7,8,12,13)]
data2<-data2[,c(1:6,8:9,7)]
write.csv(data2,file="table1.allelic.chisq.table.csv",quote=F,row.names=F)

data1<-read.table("ROI.RA3000.counts.assoc",head=T,as.is=T)
MyReferenceAllele<-subset(data1,OR<1)[,c(2,7)]
MySigSNPs<-subset(data1,P<0.05)[,c(2)]
write.table(MyReferenceAllele,file="MyReferenceAllele.txt",sep="\t",quote=F,col.names=F,row.names=F)
write.table(MySigSNPs,file="MySigSNPs.txt",sep="\t",quote=F,col.names=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)

ctabe<-function(input){
input<-as.numeric(unlist(strsplit(as.character(input),"/")))
geno.n<-input[c(3,2,1)]
geno.p<-percent(input[c(3,2,1)]/sum(input[1:3]))
dom.n<-input[1]+input[2]
dom.p<-percent((input[1]+input[2])/sum(input[1:3]))
rec.n<-input[1]
rec.p<-percent(input[1]/sum(input[1:3]))
risk.n<-2*input[1]+input[2]
risk.p<-percent((2*input[1]+input[2])/(2*sum(input[1:3])))
n<-c(geno.n,rec.n,dom.n,risk.n)
p<-c(geno.p,rec.p,dom.p,risk.p)
return(c(paste(n,"(",p,")",sep="")))
}

library("epitools")
outable<-c()
for(i in unique(data$SNP)){
temp<-subset(data,SNP==i)
mia<-temp$A1[1]
maa<-temp$A2[1]
AA<-paste(mia,mia,sep="")
AC<-paste(mia,maa,sep="")
CC<-paste(maa,maa,sep="")
rec<-paste(AA," vs ",AC,"+",CC,sep="")
dom<-paste(AA,"+",AC," vs ",CC,sep="")
A<-paste(mia,"allele",sep=" ")
TYPE<-c(CC,AC,AA,dom,rec,A)

AFF=ctabe(temp$AFF)
UNAFF=ctabe(temp$UNAFF)
data.frame(AFF,UNAFF)

REF1<-as.numeric(unlist(lapply(strsplit(as.character(AFF),split="[(]"),function(x) x[1]))[1:6])
REF2<-as.numeric(unlist(lapply(strsplit(as.character(UNAFF),split="[(]"),function(x) x[1]))[1:6])

REF<-data.frame(REF1,REF2)
OR<-"1 (reference)"
P<-1
SNP<-c()
for(j in 2:3){
xx<-as.matrix(REF[c(1,j),])
or<-round((xx[2,1]/xx[1,1])/(xx[2,2]/xx[1,2]),2)
ci95<-round(oddsratio(xx, y = NULL,method = c("fisher"), rev="col",conf.level = 0.95,correction = T,verbose = FALSE)$measure[2,2:3],2)
or<-paste(or,"(",ci95[1],"-",ci95[2],")",sep="")
p<-fisher.test(xx)$p.value
OR<-c(OR,or)
P<-c(P,p)
}
for(j in 4:5){
xx<-matrix(c(sum(REF[1:3,1])-REF[j,1],REF[j,1],sum(REF[1:3,2])-REF[j,2],REF[j,2]),2,2,byrow=F)
or<-round((xx[2,1]/xx[1,1])/(xx[2,2]/xx[1,2]),2)
ci95<-round(oddsratio(xx, y = NULL,method = c("fisher"), rev="col",conf.level = 0.95,correction = T,verbose = FALSE)$measure[2,2:3],2)
or<-paste(or,"(",ci95[1],"-",ci95[2],")",sep="")
p<-fisher.test(xx)$p.value
OR<-c(OR,or)
P<-c(P,p)
}

j=6
xx<-matrix(c(2*sum(REF[1:3,1])-REF[j,1],REF[j,1],2*sum(REF[1:3,2])-REF[j,2],REF[j,2]),2,2,byrow=F)
or<-round((xx[2,1]/xx[1,1])/(xx[2,2]/xx[1,2]),2)
ci95<-round(oddsratio(xx, y = NULL,method = c("fisher"), rev="col",conf.level = 0.95,correction = T,verbose = FALSE)$measure[2,2:3],2)
or<-paste(or,"(",ci95[1],"-",ci95[2],")",sep="")
p<-fisher.test(xx)$p.value
OR<-c(OR,or)
P<-c(P,p)
SNP<-c(SNP,rep(i,6))

outable<-rbind(outable,data.frame(TYPE,AFF,UNAFF,OR,P,SNP))
}
write.csv(outable,file="table2.model.csv",quote=F,row.names=F)
########################################################################
bedtools intersect -wao -a ~/hpc/db/hg19/refGene.hg19.V2.bed.txt -b hsa.gff3.hg19.bed | grep -v '\-1' > hsa.gff3.dbSNP153.bed
data<-read.table("ROI.RA3000.epi.cc",head=T,as.is=T)
miR<-read.table("hsa.gff3.dbSNP153.bed",head=F,as.is=T)
miR1<-miR[match(data$SNP1,miR$V9),5]
miR2<-miR[match(data$SNP2,miR$V9),5]
data<-data.frame(data,miR1,miR2)
write.csv(unique(data),file="table3.epistasis.csv",quote=F,row.names=F)






########################################################################
data<-read.table("ROI.RA3000.model",head=T,as.is=T)
GENO=na.omit(subset(data,TEST=="GENO"))
GENO=GENO[order(GENO$P,decreasing=F),]
write.csv(GENO,file="table2.GENO.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="DOM"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table3.DOM.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="REC"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table4.REC.chisq.table.csv",quote=F,row.names=F)

data<-read.table("ROI.RA3000.model",head=T,as.is=T)
input=na.omit(subset(data,TEST=="TREND"))
input=input[order(input$P,decreasing=F),]
write.csv(input,file="table6.TREND.chisq.table.csv",quote=F,row.names=F)



wget https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/cytoband.hg19.bed
perl -p -i -e 's/chr//g' cytoband.hg19.bed
bedtools intersect -wao -a hsa.gff3.hg19.bed -b cytoband.hg19.bed


