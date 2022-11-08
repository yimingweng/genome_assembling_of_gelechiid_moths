#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step3
#SBATCH -o kely_repeatmask_step3.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step3

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-e RMBlast \
-gff \
-no_is \
-lib /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely-genome-families.prefix.fa \
-dir kely_repeatmasker_step3 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step2/kely_final_assembly.fasta.masked.masked &> kely_repeatmasker_step3.out

mv kely_repeatmask_step3.out kely_repeatmasker_step3