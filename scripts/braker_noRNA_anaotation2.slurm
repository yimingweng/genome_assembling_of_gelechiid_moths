#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

genome=${1} # the masked genome
species=${2} # the species name without space (use underscore)

module load prothint/2.6.0
module load braker/2.1.6
module load busco/5.3.0
module load hmmer/3.2.1
module load diamond/2.0.9

export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# step1: run prohints to create protein database called prothint.gff
prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} /blue/kawahara/yimingweng/universal_scripts/proteins.fasta

# step2: use prothint_augustus.gff to run braker
braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} --hints=prothint_augustus.gff --softmasking --UTR=on --gff3 --cores 32 --AUGUSTUS_ab_initio

# step 3: extract the genes with full or partial support from the hints
python3 /blue/kawahara/yimingweng/universal_scripts/selectSupportedSubsets.py ./braker/augustus.hints.gtf ./braker/hintsfile.gff --fullSupport fullsupport.gff --anySupport ${species}_anysupport.gff --noSupport nosupport.gff

# step 4: redeem genes from nosupport.gff, by blasting them to the databases nr and UniProt
# step 4.1: get the amino acid sequence
/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl ${genome} nosupport.gff ${species}_nosupport.aa

# step 4.2: run diamond to blast the nosupport aa sequences to the nr and UniProt 
diamond blastp -k 1 -e 0.00001 -d /orange/kawahara/yimingweng/databases/nr.dmnd -q ${species}_nosupport.aa -f 6 -o nosupport_nr.tsv
diamond blastp -k 1 -e 0.00001 -d /orange/kawahara/yimingweng/databases/uniprot_arthropod.dmnd -q ${species}_nosupport.aa -f 6 -o nosupport_uniprot.tsv
cat nosupport_nr.tsv nosupport_uniprot.tsv | cut -d $'\t' -f 1 | sort | uniq >> transcript_list
IFS=$'\n'
for transcript in $(cat transcript_list); do cat nosupport.gff | grep -Pw ${transcript} | cut -d $'\t' -f 9 | cut -d " " -f 4 | grep -o  "[A-Za-z]_[0-9]*">> redeem_list; done

# step 4.3: add genes with blast hits to the gene model
cat ${species}_anySupport.gff | cut -d $'\t' -f 9 | cut -d " " -f 4 | grep -o  "[A-Za-z]_[0-9]*" >> redeem_list
cat redeem_list | sort | uniq >> redeem_list_uniq
sort -t "_" -k2,2 -n redeem_list_uniq >> redeem_list_uniq_sort
rm redeem_list redeem_list_uniq
IFS=$'\n'
for gene in $(cat redeem_list_uniq_sort); do cat ./braker/augustus.hints.gtf | grep -Pw ${gene} >> ${species}_braker_final.gtf; done

# step 4.4: supported_gene.gtf to generate final gene model fasta to anysupport.aa
/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl ${genome} ${species}_braker_final.gtf ${species}_braker_final.aa

# step5: run busco on the subset of the final model (remove smaller isofor)
python3 /blue/kawahara/yimingweng/universal_scripts/longest_aa.py < ${species}_braker_final.aa > ${species}_tmp.fasta

busco -f -i ${species}_tmp.fasta \
 -o ./${species}_gene_model_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32

rm ${species}_tmp.fasta ${species}_nosupport.aa