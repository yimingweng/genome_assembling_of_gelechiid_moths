#!/bin/bash

#SBATCH --job-name=Kely_minimap
#SBATCH -o Kely_minimap.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Kely_genome_project/purging/aligned.bam  \
-g /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta