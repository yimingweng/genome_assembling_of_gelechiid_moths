## **08/05/2020**
- Getting to know that we will have PacBio sequence reads for Tomato Pinworm, *Keiferia lycopersicella* and Groosefoot Groundling Moth, *Scrobipalpa atriplicella* back in a week or two. 

- The computational work will be done through the cluster machine in UF called HiperGator where an account is required.

- It might be a good idea to get familiar with this species and the genome assembly pipeline for the PacBio reads.   

<br />

## **08/08/2020**
- For the working pipeline, I can follow this paper published by [Kawahara et al., 2022](https://doi.org/10.46471/gigabyte.64). 
 This is the workflow for the paper:
   1. Genome size estimations and genome profiling ([K-Mer Counter (KMC) v.3.1.1](https://github.com/refresh-bio/KMC)) with k-mer size=21.
   2. Sequence assembly and analysis
      - Contig assembling using Hifiasm v0.13-r307
      - Assembly statistics were calculated using assembly_stats.py
      - Genome completeness was determined using BUSCO v.5.2.2 (bd10 reference Endopterygota)
      - Contamination was detected using BlobTools v1.0.
   3. Genome annotation
      - Structure Annotation: 
         - RepeatMasker (followed by RepeatModeler2) 
         - species-specific gene model training: BUSCO v.4.1.4
         - predicted genes with the homology-based gene prediction: GeMoMa v1.6.4
         - generate additional ab initio gene predictions: MAKER v3.01.03
      -Functional Annotation:
         - Add functional annotations to the predicted proteins: BlastP
         - Blast2GO for go term annotation

## **08/09/2022**
- A pilot run for the use of [K-Mer Counter (KMC) v.3.1.1](https://github.com/refresh-bio/KMC) with *Plodia interpunctella* HiFi reads.

- Get start on the Hipergator cluster (learn how to use [here](https://help.rc.ufl.edu/doc/UFRC_Help_and_Documentation))
```
ssh yimingweng@hpg.rc.ufl.edu
## type password here
```

- Find the working directory
```
# log in to the UF hpc account
[yimingweng@login6 yimingweng]$ pwd
/blue/kawahara/yimingweng
```

- Download the PacBio read file [SRR15658214](https://www.ncbi.nlm.nih.gov/sra/SRX11955122[accn]) from NCBI (the file is about 12Gb in size, and it takes about 2-3 Hrs to download using fastq-dump tool)
- For more information using slurm to submit/manage the job on the cluster mechine, check [this](https://slurm.schedmd.com/quickstart.html) out!

```
# /blue/kawahara/yimingweng

# check the latest version of sra in the machine
module spider sra
# it's sra/2.10.9
# the newer version is recently published (sra 3.0.0)

# create a slurm script
nano fastq_dump_test.slurm

#####################  script content  ####################
#!/bin/bash
#SBATCH --job-name=test_fastq_dump    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yimingweng@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=fastq_dump_test.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.9
fastq-dump --gzip  SRR15658214

###########################################################

# submit the script
sbatch fastq_dump_test.slurm

# check the job status
squeue -u yimingweng
JOBID PARTITION   NAME   USER   ST   TIME NODES NODELIST(REASON) 
44432089 hpg-milan test_fas yimingwe R  0:36  1  c0709a-s4
# once it's done, it suppose to send an email as notification

# check the result
cat fastq_dump_test.log
/blue/kawahara/yimingweng
c0702a-s1.ufhpc
Tue Aug  9 15:42:50 EDT 2022
Read 926983 spots for SRR15658214
Written 926983 spots for SRR15658214
```

- There is a concern that the adapters might remain in the CCS raw reads (see this paper by [Sim et al., 2022](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-022-08375-1)) and can potentially affect the genome assembling, although it has been claimed that the adapters are already trimmed in the PacBio Sequel system. In case the adapters are still in the reads, check the present of the adapter sequence in the fastq reads:
```
[yimingweng@login6 yimingweng]$ pwd
/blue/kawahara/yimingweng
zcat SRR15658214.fastq.gz | grep -v "@" | grep "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"  | wc -l 
# This is for the adapter, where the sequence I grep is the adapter sequence in the UniVec database of NCBI
### found 2 matches

zcat SRR15658214.fastq.gz | grep -v "@" | grep "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"  | wc -l 
# This is for the Pacific Biosciences C2 Primer
# ### found 0 match
```

## **08/10/2022**
- A pilot run for the use of [K-Mer Counter (KMC) v.3.1.1](https://github.com/refresh-bio/KMC) with *Plodia interpunctella* HiFi reads (**Continue**).
- Check the verison of KMC on the server
```
[yimingweng@login2 yimingweng]$ pwd
/blue/kawahara/yimingweng
[yimingweng@login2 yimingweng]$ module spider KMC
# $ module spider kmc/3.2.1
```

- Request a job to run KMC on the gzipped fastq file (SRR15658214.fastq.gz)
```
[yimingweng@login2 yimingweng]$ nano kmc_test.slurm

#####################  script content  ####################
#!/bin/bash
#SBATCH --job-name=test_kmc    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yimingweng@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=05:00:00               # Time limit hrs:min:sec
#SBATCH --output=test_kmc.log   # Standard output and error log
pwd; hostname; date

module load kmc/3.2.1
mkdir kmc_test
kmc -k21 SRR15658214.fastq.gz 21mers kmc_test

# Create text dump from KMC database binary format
kmc_tools transform 21mers dump 21mers.txt

# remove the temporal directory
rm -r kmc_test
###########################################################
```

- Submit the slurm job the the server
```
sbatch kmc_test.slurm
```

- I can check the job status with slurm command (i.e. squeue) or the log file
```
squeue -u yimingweng
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
44527201 hpg-defau test_kmc yimingwe  R       1:01      1 c0702a-s1

# I can monitor the progress by checking the log file with live update
tail -f test_kmc.log
/blue/kawahara/yimingweng
c0702a-s1.ufhpc
Wed Aug 10 11:49:45 EDT 2022
*************************************************
Stage 1: 100%
Stage 2: 92%
```
- Generate the histo file for plotting the result generated from KMC with the customized bash script **kmc2histo.sh**
```
bash kmc2histo.sh 21mers.txt plin_test
# the 21mers.txt is the input argument1, which should be the output from KMC

#####################  script content  ####################
#!/bin/bash
# This script is to manage the output from KMC and create the histogram that can be easily visualized through R or other programs
# This script takes 2 inputs: 1) the txt file generated by the KMC and 2) the prefix for the output file

in=${1}
out=${2}

echo "reading the input file..."
cat ${in} | cut -d $'\t' -f 2 >> out.tmp

echo "generating the temporal file..."
cat out.tmp | sort -n | uniq >> out.uniq.tmp
IFS=$'\n'
for i in $(cat out.uniq.tmp)
do
   echo "working on repeat number ${i}"
   count=$(cat out.tmp | grep -Pw ${i} | wc -c)
   echo -e "${i}\t${count}" >> ${out}.histo
done
rm *tmp
```
- Transfer the histo file to the local machine, so that I can use R to visualize the result
```
# yiming@DESKTOP-H41N7NT:/mnt/c/Users/wengz/Dropbox/postdoc/Kely_genome_project$
scp -P 22 yimingweng@hpg.rc.ufl.edu:/home/yimingweng/plin.histo ./
```


- plot the histogram using R script **kmer_plot.r**
```
######## R environment ########
path <- getwd()

myhisto <- read.table(paste0(path, "/plin.histo"), header=F, sep="\t")
x11() # this only works for windows R I guess
plot(myhisto,
     type="l",
     xlab="kmer count",
     ylab="frequency") # plots the line from the data points 

# replot the histogram
highest <- myhisto[which(myhisto$V2 == max(myhisto[10:200,2])),]

x11() # this only works for windows R I guess
plot(myhisto[10:200,],
     type="l",
     xlab="kmer count",
     ylab="frequency") # plots the line from the data points 
points(myhisto[10:200,]) # plot the data points
text(highest[,1], highest[,2], highest$V1,
     cex=1, pos=2,col="black")
```