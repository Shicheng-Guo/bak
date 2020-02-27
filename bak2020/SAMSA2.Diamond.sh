#####################################################
## diamond and SAMSA2
#####################################################
mkdir ~/hpc/tools/diamond
cd ~/hpc/tools/diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.30/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr

diamond makedb --in nr.faa -d nr
diamond blastx -d nr -q reads.fna -o matches.m8

git clone https://github.com/transcript/samsa2.git
cd setup_and_test
sh ./package_installation.bash
