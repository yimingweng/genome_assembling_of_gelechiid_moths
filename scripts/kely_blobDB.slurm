#!/bin/bash

#SBATCH --job-name=Kely_blobDB
#SBATCH -o Kely_blobDB.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 8:00:00
#SBATCH -c 8

module load blobtools/2.2

# create blobDB
blobtools create \
--fasta /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_15X_120X.fasta \
--cov /blue/kawahara/yimingweng/Kely_genome_project/blobplot/purging_aligned.bam \
Kely_purge_DB

# get the node.dmp file from NCBI, and put it in ./Kely_purge_DB
cd ./Kely_purge_DB
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xzf new_taxdump.tar.gz
cd ..

# add the tax and blast results
blobtools add \
--taxdump /blue/kawahara/yimingweng/Kely_genome_project/blobplot/Kely_purge_DB/ \
--taxid 1511202 \
--hits /blue/kawahara/yimingweng/Kely_genome_project/blobplot/kely_purge_assembly.nt.mts1.hsp1.1e25.megablast.out \
Kely_purge_DB

# create summary tsv file
blobtools filter --table Kely_purge_summary.tsv /blue/kawahara/yimingweng/Kely_genome_project/blobplot/Kely_purge_DB/