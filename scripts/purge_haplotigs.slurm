#!/bin/bash

#SBATCH --job-name=Kely_purge_haplotigs
#SBATCH -o Kely_purge_haplotigs.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta  \
-c /blue/kawahara/yimingweng/Kely_genome_project/purging/coverage_stats.csv \
-o Kely_purge_5X_120X