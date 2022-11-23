#!/bin/bash
#SBATCH --job-name=%x_braker_prothint_%j
#SBATCH --output=%x_braker_prothint_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1} # the masked genome
species=${2} # the species name without space (use underscore)

module load prothint/2.6.0
module load braker/2.1.6
module load busco/5.3.0
module load hmmer/3.2.1

export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# step1: run prohints to create protein database called prothint.gff
prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} /blue/kawahara/yimingweng/universal_scripts/proteins.fasta

# step2: use prothint_augustus.gff  to run braker
braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} --hints=prothint_augustus.gff --softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio

# step 3: extract the genes with full or partial support from the hints
python3 /blue/kawahara/yimingweng/universal_scripts/selectSupportedSubsets.py ./braker/augustus.hints.gtf ./braker/hintsfile.gff --fullSupport fullsupport.gff --anySupport ${species}_prohint_braker_final.gff --noSupport nosupport.gff

# step 4: supported_gene.gtf to generate final gene model fasta called "braker_supported_gene.aa"
/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl ${genome} ${species}_prohint_braker_final.gff ${species}_anysupport_aa.fa

# step 5: run busco on the subset of the final model (remove smaller isofor)
python3 /blue/kawahara/yimingweng/universal_scripts/longest_aa.py < ${species}_anysupport_aa.fa > ${species}_tmp.fasta

busco -f -i ${species}_tmp.fasta \
 -o ./${species}_gene_model_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32

rm ${species}_tmp.fasta

# step 6: get the final gft 
cat ${species}_prohint_braker_final.gff | cut -d $'\t' -f 9 | cut -d " " -f 2 | cut -d "\"" -f 2 | grep -v "#" | sort -t "g" -k2,2 -n | uniq >> transcript_list
touch ${species}_prohint_braker_final.gtf
IFS=$'\n'
for transcript in $(cat transcript_list)
do
    gene=$(echo ${transcript} | cut -d "." -f 1)
    gene_line=$(cat ./braker/augustus.hints.gtf | grep -Pw "gene" | grep -Pw "${gene}")
    match=$(cat ${species}_prohint_braker_final.gtf | grep -Pw "gene" | grep -Pw "${gene}")
    if [[ -z ${match} ]]
    then
        echo ${gene_line} >> ${species}_prohint_braker_final.gtf
    else
        echo "gene is already in the gft, skip"
    fi
    cat ./braker/augustus.hints.gtf | grep -Pw "${transcript}" >> ${species}_prohint_braker_final.gtf
done