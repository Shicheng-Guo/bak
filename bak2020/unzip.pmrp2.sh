wget https://imputationserver.sph.umich.edu/share/results/46ef612c40ec44201cc7f9e9ae4a005e/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/96b629e0627d48c936b175bdf75004b2/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/f7a64d692e74e6f45d3391f30536fce9/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/c64a8035e84d05a7190b6c77fe5b823f/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/6a560e42e07ae33b4d7621568dc37ae4/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/c7a6fc8599f2509a9b99266727c74d73/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/5148e156887471878650d402b1bf4cbf/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/1bba3f4004b90fb6d5b169e572fe449c/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/d148dbd3b6dce42a1c3afdabc2098a6/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/598d45e4d38fce1b31acad817d67f39e/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/ef0d728e8ceb136a8d16998b7a19e2a8/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/52dd42d0b07e922f087d92bbb77a5c8a/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/ec77fb9b950dbb9dd2e1493dd53bef74/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/ddafbd95610e1093aeb3ae50d0582e05/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/fe4faeef120657ac1170c0b6b6bfccce/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/1b888a9e2e80ab9fe12f391087b690ed/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/144e4596526eb6f807595eb64c7eb53e/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/1c2358cb9055c296171fe946f4dc8ceb/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/59700561aa23dc9a4894ca7bc37368f5/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/6ffaa7b7c35e3cd820ac96198457bd58/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/c5924f254f0437959724975c37b040a6/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/ff3aa625bab98c0b574d3405f37d0b11/chr_9.zip
wget https://imputationserver.sph.umich.edu/share/results/fcbc5c38e5104812d8bd4a82d3643fe5/chr_X.zip

mkdir temp
for i in {1..23} X
do
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=4 >> $i.job
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
echo unzip -P \'L9n6EY\:nNsxyFI\' chr_$i.zip  >> $i.job
qsub $i.job
done

