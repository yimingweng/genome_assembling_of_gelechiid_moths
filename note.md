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

- get start on the Hipergator cluster
```

```

- find the working directory
```
# log in to the UF hpc account
[yimingweng@login6 yimingweng]$ pwd
/blue/kawahara/yimingweng
```

- download the PacBio read file [SRR15658214](https://www.ncbi.nlm.nih.gov/sra/SRX11955122[accn]) from NCBI
- for more information using slurm to submit/manage the job on the cluster mechine, check [this](https://slurm.schedmd.com/quickstart.html) out!

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
JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)                                        44432089 hpg-milan test_fas yimingwe  R       0:36      1 c0709a-s4
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
# ### found o match
```