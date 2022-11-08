#!/bin/bash

#SBATCH --job-name=kely_megablast_nt
#SBATCH -o kely_megablast_nt.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 72:00:00
#SBATCH -c 64

module load ncbi_blast/2.10.1

blastn \
-task megablast \
-query /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_15X_120X.fasta \
-db nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 64 \
-evalue 1e-25 \
-out kely_purge_assembly.nt.mts1.hsp1.1e25.megablast.out