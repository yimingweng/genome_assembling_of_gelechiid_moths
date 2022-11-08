#!/bin/bash

#SBATCH --job-name=Kely_minimap_purge
#SBATCH -o Kely_minimap_purge.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 8:00:00
#SBATCH -c 8

module load minimap/2.21
module load samtools/1.15

minimap2 -t 8 -ax map-pb /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_15X_120X.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o purging_aligned.bam -T tmp.ali