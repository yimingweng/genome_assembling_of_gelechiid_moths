#!/bin/bash

#SBATCH --job-name=%x_interproscan_%j
#SBATCH -o %x_interproscan_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 120:00:00
#SBATCH -c 32

module load iprscan/5.57

protein=${1}
prefix=${2}

cat ${protein} | sed 's/\*//g' > ${prefix}.fasta
interproscan.sh -i ${prefix}.fasta -cpu 32 -f tsv -goterms