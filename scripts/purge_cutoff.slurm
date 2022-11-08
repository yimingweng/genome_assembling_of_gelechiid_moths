#!/bin/bash

#SBATCH --job-name=purge_cutoff
#SBATCH -o purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Kely_genome_project/purging/aligned.bam.gencov  \
-l 5  \
-m 33  \
-h 120  \
-o coverage_stats.csv \
-j 80 \
-s 80