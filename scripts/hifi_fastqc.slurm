#!/bin/bash
#SBATCH --job-name=%x_fastqc_%j
#SBATCH --output=%x_fastqc_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=02:00:00

# sometimes the issue with perl version occur, then you see the error "symbol lookup error"
# If so, try different version of perl and run fastqc

## load order version of perl
module load perl/5.20.0

## load fastqc module
module load fastqc/0.11.7

input=${1}

## run fastqc on the data
fastqc ${input}