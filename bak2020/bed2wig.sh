#!/bin/bash
cp /mnt/gluster/ENCFF003JVR.bed.bedgraph ./
cp /mnt/gluster/hg38.chrom.sizes ./
wigToBigWig ENCFF003JVR.bed.bedgraph hg38.chrom.sizes ENCFF003JVR.bw
cp ENCFF003JVR.bw /mnt/gluster/nu_guos
rm ENCFF003JVR.bed.bedgraph
rm ENCFF003JVR.bw
rm hg38.chrom.sizes
