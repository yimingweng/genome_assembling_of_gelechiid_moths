#!/bin/bash
#SBATCH --job-name=%x_gfacs_default
#SBATCH --output=%x_gfacs_default.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load perl/5.24.1

format=${1}
model=${2}
genome=${3}
prefix=${4}
outdir=${5}

perl /blue/kawahara/yimingweng/universal_scripts/gFACs/gFACs.pl \
-f ${format} \
-p ${prefix} \
    --statistics \
    --statistics-at-every-step \
--distributions \
    exon_lengths 10 \
    intron_lengths 15 \
    CDS_lengths 20 \
    gene_lengths 100 \
    exon_position \
    exon_position_data \
    intron_position \
    intron_position_data \
--fasta ${genome} \
    --nt-content \
--create-gff3 \
--get-protein-fasta \
--create-gtf \
-O ${outdir} \
${model}