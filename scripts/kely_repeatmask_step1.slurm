#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step1
#SBATCH -o kely_repeatmask_step1.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step1

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-noint \
-no_is \
-dir kely_repeatmasker_step1 \
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta &> kely_repeatmasker_step1.out