#!/bin/bash
#SBATCH -J call_peak
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=50000
#SBATCH -e errput_%J
#SBATCH -o output_%J

source activate archr_2
file=$1
sample_genome=/cluster/home/huangyiming/phast/fa_size/mm10.chrom.sizes
sort -k 1,1 -k2,2n ${file}/${file}_blacklistrm.bed > ${file}/${file}_blacklistrm.sort.bed
bedtools genomecov -i ${file}/${file}_blacklistrm.sort.bed -bg -g ${sample_genome} > ${file}/${file}.bedGraph
bedGraphToBigWig ${file}/${file}.bedGraph ${sample_genome} ${file}/${file}.bw