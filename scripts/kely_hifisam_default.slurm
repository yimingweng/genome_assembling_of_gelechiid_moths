#!/bin/bash
#SBATCH --job-name=hifiasm_reassemble_kely
#SBATCH -o hifiasm_reassemble_kely.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH -c 32
#SBATCH --mem-per-cpu=10gb
#SBATCH -t 30:00:00
#SBATCH --account=kawahara
#SBATCH --qos=kawahara-b

module load ufrc
module load hifiasm

# to optimize the power of purging, try "-l 3"
hifiasm -o /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam.asm \
        -l 3 \
        -t 32 \
        /blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz

# convert the output file to fasta
awk '/^S/{print ">"$2;print $3}' /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam.asm.bp.p_ctg.gfa > /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta
