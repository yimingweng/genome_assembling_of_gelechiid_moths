#!/bin/bash

#SBATCH --job-name=kely_repeatmodeler2.slurm
#SBATCH -o kely_repeatmodeler2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 08:00:00
#SBATCH -c 30

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmodeler/2.0
module load seqkit/2.0.0

genome=${1} # full path to the genome assembly
outprefix=${2}

# build the RM2 database for the genome
BuildDatabase -name kely_genome ${genome}

# run RepeatModeler with the database
RepeatModeler -database kely_genome -pa 10 -LTRStruct >& ${outprefix}.out

cat kely_genome-families.fa | seqkit fx2tab | awk '{ print "Kely_1.0_"$0 }' | seqkit tab2fx > kely_genome-families.prefix.fa