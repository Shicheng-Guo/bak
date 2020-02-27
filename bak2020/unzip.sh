

wget https://imputationserver.sph.umich.edu/share/results/65d218b7982673096c80c100f823e510/chr_1.zip &
wget https://imputationserver.sph.umich.edu/share/results/aa1fa3b21c7d2a99e44ed2e9e58dc45e/chr_10.zip  &
wget https://imputationserver.sph.umich.edu/share/results/4c8e856dc68ce9c1f5a27453c0bd09f6/chr_11.zip &
wget https://imputationserver.sph.umich.edu/share/results/bb6d350c35f0b66032c10b3163894d25/chr_12.zip &
wget https://imputationserver.sph.umich.edu/share/results/a497325321bf82458897ac4386b6a3ba/chr_13.zip &
wget https://imputationserver.sph.umich.edu/share/results/952ebe6516a7f820862942843da21e08/chr_14.zip &
wget https://imputationserver.sph.umich.edu/share/results/b27065f7d08208f05b8b836ceb0c7b30/chr_15.zip &
wget https://imputationserver.sph.umich.edu/share/results/4af52101004ce7ba2bc53d9710516e5e/chr_16.zip &
wget https://imputationserver.sph.umich.edu/share/results/3ba2e5d1b186b999618df8c5195f67d4/chr_17.zip &
wget https://imputationserver.sph.umich.edu/share/results/47e51e23b1ec8b3a42ee0fae0c3eddf9/chr_18.zip &
wget https://imputationserver.sph.umich.edu/share/results/79a32ec4e41032bdadd9c71ed9e01026/chr_19.zip &
wget https://imputationserver.sph.umich.edu/share/results/388f3945331eefe98bcdbe64242a4c6/chr_2.zip &
wget https://imputationserver.sph.umich.edu/share/results/8b9b0517faa5fdd8796ee174545f75ea/chr_20.zip &
wget https://imputationserver.sph.umich.edu/share/results/290c3c5364b01bd8bc4ca2b58ac38067/chr_21.zip &
wget https://imputationserver.sph.umich.edu/share/results/90d4670ced179a6e594b65e13c6c972a/chr_22.zip &
wget https://imputationserver.sph.umich.edu/share/results/f2ca251cc17057a34832fb9c6b01385b/chr_3.zip &
wget https://imputationserver.sph.umich.edu/share/results/c17f5e206e5e20e6130c10a77a18e63f/chr_4.zip &
wget https://imputationserver.sph.umich.edu/share/results/65f3291df7ee40d6ef12a5cee27851b/chr_5.zip &
wget https://imputationserver.sph.umich.edu/share/results/54c42835635bcdc8b8f0634e51e8fcf1/chr_6.zip &
wget https://imputationserver.sph.umich.edu/share/results/d1509e1d32a2e1f73e901d59bfed71c5/chr_7.zip &
wget https://imputationserver.sph.umich.edu/share/results/92c07addbc9da6e89390e8dd7b782b5/chr_8.zip &
wget https://imputationserver.sph.umich.edu/share/results/9ba611afa53b403df13ab963715111a2/chr_9.zip &
wget https://imputationserver.sph.umich.edu/share/results/785e3e9ef0ec5df953c7f03ca216681c/chr_X.zip &

mkdir temp
for i in {1..23} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo plink --bfile $plink --chr $i --recode vcf-iid --out $plink.chr$i >> $i.job
# echo java -jar ./conform-gt.24May16.cee.jar gt=$plink.chr$i.vcf chrom=$i ref=~/hpc/db/hg19/beagle/EUR/chr$i.1kg.phase3.v5a.dedup.norm.EUR.vcf.gz out=$plink.chr$i.beagle >> $i.job
# echo bcftools norm -m-any chr$i.1kg.phase3.v5a.EUR.vcf.gz -Oz -o chr$i.1kg.phase3.v5a.dedup.EUR.vcf.gz  >>$i.job
# echo bcftools norm -d all chr$i.1kg.phase3.v5a.dedup.EUR.vcf.gz -Oz -o chr$i.1kg.phase3.v5a.dedup.norm.EUR.vcf.gz  >>$i.job
# echo plink --bfile S_Hebbring_Unr.Guo --recode vcf --chr $i --snps-only just-acgt --out ./beagle/S_Hebbring_Unr.Guo.Forward.chr$i >> chr$i.job
echo unzip -P \'T8gBwPI8\$dUilF\' chr_$i.zip  >> $i.job
qsub $i.job
done