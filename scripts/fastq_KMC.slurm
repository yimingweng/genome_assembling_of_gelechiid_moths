#!/bin/bash
#SBATCH --job-name=%x_kmc_%j
#SBATCH --output=%x_kmc_%j.log
#SBATCH --mail-type=END,FAI
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=05:00:00 

input=${1}
prefix=${2}
kmer=${3}

module load kmc/3.2.1
mkdir kmc_tmp

kmc -k${kmer} ${input} ${prefix}_${kmer}mers kmc_tmp

# Create text dump from KMC database binary format
kmc_tools transform ${prefix}_${kmer}mers dump ${prefix}_${kmer}mers.txt

# Create histogram from the txt file
kmc_tools transform ${prefix}_${kmer}mers histogram ${prefix}_${kmer}mers.histo

# remove the temporal directory
rm -r kmc_tmp