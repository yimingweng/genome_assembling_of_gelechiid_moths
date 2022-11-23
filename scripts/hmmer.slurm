#!/bin/bash

#SBATCH --job-name=%x_hmmer_%j
#SBATCH -o %x_hmmer_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 32

module load hmmer/3.2.1

seqfile=${1}
species=${2}

hmmscan --tblout ${species}_hmmer /orange/kawahara/yimingweng/databases/Pfam-A.hmm ${seqfile}
[yimingweng@login1 hmmer]$ cp hmmer.slurm /blue/kawahara/yimingweng/universal_scripts/