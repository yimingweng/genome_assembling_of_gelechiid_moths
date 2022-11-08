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
- For more information using slurm to submit/manage the job on the cluster machine, check [this](https://slurm.schedmd.com/quickstart.html) out!

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
- Check the version of KMC on the server
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
minimap2
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
###########################################################
```
- Transfer the histo file to the local machine, so that I can use R to visualize the result
```
# yiming@DESKTOP-H41N7NT:/mnt/c/Users/wengz/Dropbox/postdoc/Kely_genome_project$
scp -P 22 yimingweng@hpg.rc.ufl.edu:/home/yimingweng/plin.histo ./
```


- Plot the histogram using my R script **kmer_plot.r**
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

- Here is this plot:
<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/kmer_plots/plin_kmer_density_plot.jpg?raw=true?raw=true">

<br />

## **08/31/2022**  
**\# PacBio read handling**  
**\# genome assembly**  

The first genome sequence reads are now available. So let's play with this giant package (2.1T).
1. Under the path I desired to store the reads (usually in `orange/`), run the command provided by the sequencing company or UF-ICBR
```
[yimingweng@login5 2_B03]$ pwd
/orange/kawahara/yimingweng/raw_read_storage/NS2846/2_B03
# run the command provided by the sequencing facility
# this takes hours to finish, so make sure the ssh can stay linked for this long period of time to allow the data transition finish
/orange/icbrdd/receive NS2846 

# when the data transition is done, a folder shows up and there are lots of bam files in there
# the one with extension of "subreads.bam" is the one with clean reads 
[yimingweng@login5 2_B03]$ ls
hifi_reads				                       m64219e_220809_225049.subreadset.xml
m64219e_220809_225049.baz2bam_1.log         m64219e_220809_225049.subreadset.xml.bak
m64219e_220809_225049.scraps.bam	           m64219e_220809_225049.transferdone
m64219e_220809_225049.scraps.bam.pbi	     subreadsTohifi-2.sh
m64219e_220809_225049.sts.xml	      	     subreadsTohifi-bi.sh
m64219e_220809_225049.subreads.bam	        subreadsTohifi.sh
m64219e_220809_225049.subreads.bam.pbi	     tmp-file-ce9f726a-c3c4-46e6-a88f-64c178b33ed6.txt
```
2. Although there is a fastq file in the `hifi_reads` folder, It is worth to keep the record of the filtering steps. So here is the filtering steps when converting the files file from bam to fastq
```
#####################  script content  ####################
bamtools filter -in m64219e_220809_225049.subreads.bam -out hifireads/hifi_reads.bam -tag "rq":>=0.99

bamtools filter -in m64219e_220809_225049.subreads.bam -out hifireads-1/hifi_reads.bam -tag "rq":>=0.99

ccs --hifi-kinetics --min-rq 0.99 --report-file hifireads-2/hifi_report.txt m64219e_220809_225049.subreads.bam hifireads-2/hifi_reads.bam
###########################################################
```

3. To further check the read quality, use FASTQC to get a look at the quality summary.
```
[yimingweng@login2 raw_reads]$ sbatch hifi_fastqc.slurm Keiferia_lycopersicella_ccs.fastq.gz

#####################  script content  ####################
#!/bin/bash
#SBATCH --job-name=hifi_read_fastqc
#SBATCH --output=hifi_read_fastqc.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=02:00:00

# sometimes the issue with perl version occur, then you see the error "symbol lookup error"
# If so, try different version of perl and run fastqc

## load order version of perl
module load perl/5.20.0

## load fastqc module
module load fastqc/0.11.7

input=${1}

## run fastqc on the data
fastqc ${input}
###########################################################
```

- It seems ok with the data, as the GC content distribution is a bit sharp. It is unlikely to be the adaptor contamination as they should have been trimmed and also there is nothing detected in the overrepresented sequences.

- collect results and move them into: `/blue/kawahara/yimingweng/Kely_genome_project/raw_reads/fastqc`
```
[yimingweng@login2 raw_reads]$ mkdir fastqc
[yimingweng@login2 raw_reads]$ mv *fastqc* fastqc/
```

4. Check the basic statistics of the raw read data.
```
[yimingweng@login2 raw_reads]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/raw_reads

zcat Keiferia_lycopersicella_ccs.fastq.gz | grep "/ccs$" | wc -l
3839580 # there are 3,839,580 reads (~3.8 million reads) in the files

zcat Keiferia_lycopersicella_ccs.fastq.gz | grep -v [^ATCG] | wc -c
21133733487 # there are 21133733487 bps in the data (~ 21X if genome size is 500mbp, ~35X if genome is 300mbp)
```

5. Run KMC with the fastq read file
```
[yimingweng@login2 genome_size]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/genome_size

sbatch fastq_KMC.slurm /blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz Kely 21

#####################  script content  ####################
#!/bin/bash
#SBATCH --job-name=Kely_kmc
#SBATCH --output=Kely_kmc.log
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
###########################################################
```

6. Get the kmer count distribution using [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/). The result can be found [here](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=mYfpkD8wYxxcIy48XpO0)
- The estimated genome size is about 302 Mbp
- The heterozygosity is fairly low, about 1%
- The mean coverage is about 24.7X, slightly lower than expected
- This is possibly due to the repeat regions at the 100X peak which is odd, as it seems that part of the genome is **tetraploid**, as the repeats are ~25/~50/~100.
![](http://qb.cshl.edu/genomescope/genomescope2.0/user_data/mYfpkD8wYxxcIy48XpO0/linear_plot.png)
![](http://qb.cshl.edu/genomescope/genomescope2.0/user_data/mYfpkD8wYxxcIy48XpO0/transformed_linear_plot.png)

7. Let's try K=27 and see if the 100X peak is retained.
```
[yimingweng@login2 genome_size]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/genome_size

sbatch fastq_KMC.slurm /blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz Kely 27
```
The result of K=27 can be found [here](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=3w5G4b8elEYO96CcJdTX)
- The estimated genome size is about 318 Mbp
- The heterozygosity is still low, about 0.9%
- The mean coverage is about 24.5X
- The 100X peak is still there, which means that the repeats have some pattern or some duplication occurred. This is not novel, check this [example by Nagy et al](https://link.springer.com/article/10.1186/s12864-021-07627-w).
![](http://qb.cshl.edu/genomescope/genomescope2.0/user_data/3w5G4b8elEYO96CcJdTX/linear_plot.png)
![](http://qb.cshl.edu/genomescope/genomescope2.0/user_data/3w5G4b8elEYO96CcJdTX/transformed_linear_plot.png)

8. Assemble the genome with hifisam
```
[yimingweng@login6 assemblies]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies

mkdir kely_hifisam_default
sbatch kely_hifisam_default.slurm

###########################  script content  ###########################
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
########################################################################
```
<span style="color:red">NOTE: the output files of KMC are pretty large, so delete them and just keep the histogram files. </span>

<br />

## **09/01/2022**  
**\# BUSCO**  
**\# purge_haplotigs pipeline**  
**\# blobplot**  

The genome assembly is done and let's look at the assembly statistics.

1. Use the previously developed python script "`assemblystats.py`" to get the assembly statistics.
```
[yimingweng@login5 kely_hifisam_default]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default

# allow the script to be excuted
[yimingweng@login5 kely_hifisam_default]$ chmod +x assemblystats.py

# run the script with the input of the assembly
./assemblystats.py Kely_hifisam.asm.bp.p_ctg.gfa Kely_hifisam_default_stats

# keep  the script for future reference and use
[yimingweng@login5 kely_hifisam_default]$ cp ./assemblystats.py /blue/kawahara/yimingweng/scripts/

#
[yimingweng@login5 kely_hifisam_default]$ cat Kely_hifisam_default_stats
{
  "Contig Stats": {
    "L10": 1,
    "L20": 3,
    "L30": 6,
    "L40": 8,
    "L50": 12,
    "N10": 24805861,
    "N20": 22264078,
    "N30": 16479331,
    "N40": 15732994,
    "N50": 14121271,
    "gc_content": 38.88312173637028,
    "longest": 29025318,
    "mean": 3928602.3017241377,
    "median": 460798.0,
    "sequence_count": 116,
    "shortest": 6552,
    "total_bps": 455717867
  },
  "Scaffold Stats": {
    "L10": 1,
    "L20": 3,
    "L30": 6,
    "L40": 8,
    "L50": 12,
    "N10": 24805861,
    "N20": 22264078,
    "N30": 16479331,
    "N40": 15732994,
    "N50": 14121271,
    "gc_content": 38.88312173637028,
    "longest": 29025318,
    "mean": 3928602.3017241377,
    "median": 460798.0,
    "sequence_count": 116,
    "shortest": 6552,
    "total_bps": 455717867
  }
}

###########################  script content  ###########################
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
########################################################################
```
2. The result looks pretty good. It has **116 contigs** with N50 of **14,121,271 bps**. I have a little concern about the assembled size which is about **455 Mbps**, larger than the estimate by GenomeScope (318 Mbp). This could be owing to the presense of duplication regions (the third peak in the kmer count distribution)? Maybe BUSCO can help clarify the situation. If the duplication rate is high, maybe I should try to run the purge haplotig pipeline.
- run busco for this genome
```
[yimingweng@login5 kely_hifisam_default]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default

# copy the Augustus contig folder to here
[yimingweng@login5 kely_hifisam_default]$ cp -r /blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/ ./

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_busco
#SBATCH -o Kely_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 3:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0

# run busco command
busco -f -i /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta \
 -o ./busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m genome -c 6 --augustus
########################################################################

# check the results
[yimingweng@login2 busco_out]$ cat short_summary.specific.endopterygota_odb10.busco_out.
cat: short_summary.specific.endopterygota_odb10.busco_out.: No such file or directory
        C:98.1%[S:95.8%,D:2.3%],F:0.5%,M:1.4%,n:2124
        2082    Complete BUSCOs (C)
        2034    Complete and single-copy BUSCOs (S)
        48      Complete and duplicated BUSCOs (D)
        10      Fragmented BUSCOs (F)
        32      Missing BUSCOs (M)
        2124    Total BUSCO groups searched
```
3. The result of busco looks good, it has the completeness (C) of **98.1%** but like it is expected, the duplication rate is a bit high (it's **2.3%** and usually I see duplication rate between 0.5%-2%). So let's go over the [purge_haplotigs pipeline](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/). There are several steps in this pipeline, let's do it step by step.
-  step 1: map the raw reads (subreads) to the target genome assembly.

```
[yimingweng@login2 kely_hifisam_default]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default

# prepare the input files
[yimingweng@login2 kely_hifisam_default]$ gzip -c Kely_hifisam_default.fasta > Kely_hifisam_default.fasta.gz
```


```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_minimap
#SBATCH -o Kely_minimap.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Kely_genome_project/purging/aligned.bam  \
-g /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta  
########################################################################
```
- step 2: determine the cutoffs of coverage for the contigs based on the coverage histogram. The cutoffs I would like to apply is min=5 and max=120. I will use this values to purge the contigs.
<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/aligned.bam.histogram.png?raw=true?raw=true">

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=purge_cutoff
#SBATCH -o purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Kely_genome_project/purging/aligned.bam.gencov  \
-l 5  \
-m 33  \
-h 120  \
-o coverage_stats.csv \
-j 80 \
-s 80
########################################################################
```

- step 3: Run the purging pipeline.

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_purge_haplotigs
#SBATCH -o Kely_purge_haplotigs.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 3:00:00
#SBATCH -c 4

module load minimap/2.21
module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta  \
-c /blue/kawahara/yimingweng/Kely_genome_project/purging/coverage_stats.csv \
-o Kely_purge_5X_120X
########################################################################
```

- The resulting genome assembly has only **62 contigs** (originally it has 116 contigs), but the assembled size doesn't reduce much (**455 Mbps >> 444 Mbps**). It seems that it has removed most of the smallest contigs even though the cutoffs were not about size but coverages.

4. Rerun BUSCO on the purged genome.
```
[yimingweng@login6 Kely_purge_5X_120X]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/busco/Kely_purge_5X_120X

sbatch universal_run_busco.slurm /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta Kely_purge_5X_120X

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=universal_run_busco
#SBATCH -o universal_run_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 3:00:00
#SBATCH -c 12

fasta=${1} # the full path of genome fasta file
outprefix=${2} # the name desired to be for the putput

# copy the contig files from `/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/
cp -r /blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus ./
cp -r /blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/config.ini ./

export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

#load busco, make sure this is the latest version
module load busco/5.3.0

mkdir ${outprefix}

#run busco command
busco -f -i ${fasta} \
 -o ./${outprefix} \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m genome -c 6 --augustus

rm -rf Augustus config.ini
########################################################################
```
- The result shows that the genome has certain improvement. Now the BUSCO is:
BUSCO: C:98.0%[S:96.9%,D:1.1%],F:0.5%,M:1.5%,n:2124. However, it will be better to further look at the GC content to identify the non-target sequences. So let's do blobplot.


5. Using [blobplot](https://blobtools.readme.io/docs) to identify the non-target sequences in the genome assembly. See the instruction [here](https://blobtools.readme.io/docs/my-first-blobplot) to do the work.
- create the database (blastdb) for [hits](https://blobtools.readme.io/docs/taxonomy-file) file, one of the required input to make blobplot. <span style="color:red"> After testing `blastn` function without the downloaded database, the local nt database is not required. So I have removed the database as it takes too much storing space <sapn />
```
[yimingweng@login2 blastdb]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot/blastdb

# download the perl script (update_blastdb.pl) from NCBI
# make it excusible
[yimingweng@login2 blastdb]$ chmod +x  update_blastdb.pl

# to make it work, make sure you load the perl module
module load perl/5.24.1

#test this perl script
./update_blastdb.pl --help

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=ncbi_nt_download
#SBATCH -o ncbi_nt_download.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=1gb
#SBATCH -t 12:00:00
#SBATCH -c 24

chmod +x  update_blastdb.pl
module load perl/5.24.1

perl ./update_blastdb.pl --passive --decompress nr
########################################################################
```
2. Although it is not sure I will have to download the database, I can take the next step to blast the contigs to the nt database. 

```
sbatch kely_megablast_nt.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=kely_megablast_nt
#SBATCH -o kely_megablast_nt.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 72:00:00
#SBATCH -c 64

module load ncbi_blast/2.10.1

blastn \
-task megablast \
-query /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta \
-db nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 64 \
-evalue 1e-25 \
-out kely_purge_assembly.nt.mts1.hsp1.1e25.megablast.out
########################################################################
```

<br />

## **09/03/2022**  
**\# blobplot**  
1. Now prepare another input files for blobplot: the mapped reads file. Note that I have used minimap2 to do the similar work but that was on the original assembly. So I will have to repeat this work on the purged assembly again for the blobplot.

```
[yimingweng@login5 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

sbatch kely_minimap_purge.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=kely_minimap_purge
#SBATCH -o kely_minimap_purge.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 8:00:00
#SBATCH -c 8

module load minimap/2.21
module load samtools/1.15

minimap2 -t 8 -ax map-pb /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o purging_aligned.bam -T tmp.ali
########################################################################
```

## **09/05/2022**  
**\# blobplot** 
1. With the three input files ready (the genome assembly, the reads mapped to the assembly in bam format, and the hit file by blasting the contig to the NCBI nt database), let's run blobplot create function to creating a blobDB. Note that the example command in the [instruction](https://blobtools.readme.io/docs/my-first-blobplot) at this step is not very useful (at least to me), so consider googling around if encountering errors.

```
[yimingweng@login6 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

sbatch kely_blobDB.slurm

###########################  script content  ###########################
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
--fasta /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta \
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
########################################################################
```
2. Once the blobtool works without error, we should see an output file generated in the current working directory (not in the database folder). Here it's called Kely_purge_summary.tsv
- use local computer (my dell laptop) to draw the plot with R (I am drawing the blobplot by myself because I can't find the command line function working in blobtools v2 to draw the plot, and it looks like an [unfinished work](https://github.com/blobtoolkit/blobtoolkit/issues/16). Or it seems reply on the GUI tool called BlobToolKit Viewer).

```
yiming@DESKTOP-H41N7NT:/mnt/c/Users/wengz/Dropbox/postdoc/Kely_genome_project$ pwd
/mnt/c/Users/wengz/Dropbox/postdoc/Kely_genome_project

scp  yimingweng@hpg.rc.ufl.edu:/blue/kawahara/yimingweng/Kely_genome_project/blobplot/Kely_purge_summary.tsv ./
```

- Run the Rscipt to generate the plot


```
###### R environment ######
# this script is used to plot the blobplot for Keiferia genome
# please modify this script for the use on the other cases
library(ggplot2)


path <- getwd()
contig_dat <- read.table(paste0(path, "/Kely_purge_summary.tsv"), header=T, sep="\t")

x11()
ggplot(contig_dat, aes(x=gc, y=purging_aligned_cov, size=length, col=bestsumorder_phylum)) + 
  scale_color_manual(values=c("blue", "darkgreen")) +
  geom_point(alpha=0.2) +
  scale_size(range = c(2, 10), name="length") +
  scale_x_continuous(name="Speed of cars", limits=c(0, 1)) +
  ylab("Contig coverage") +
  xlab("GC proportion") +
  guides(size = "none") +
  guides(color=guide_legend(override.aes = list(size=3), title="best blast phylum")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```
    
- Here is the result, there is a very small contig likely to be from a plant. I think it needs to be removed from the genome assembly.
<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/blobplot/Kely_purge_blobplot.jpeg?raw=true?raw=true">

## **09/06/2022**  
**\# blobplot**

1. It would be interesting to look at the blobplot for the original (default) genome assembly too. So repeat the work for the original assembly to see what were trimmed by the purging step for future reference.
- blast the original assembly

```
[yimingweng@login5 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplo

sbatch megablast_nt.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta kely_hifisam_default

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=kely_megablast_nt
#SBATCH -o kely_megablast_nt.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 72:00:00
#SBATCH -c 64

module load ncbi_blast/2.10.1

assembly=${1}
outprefix=${2}

blastn \
-task megablast \
-query ${assembly} \
-db nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 64 \
-evalue 1e-25 \
-out ${outprefix}.megablast.out
########################################################################
```

- make the plot database

```
[yimingweng@login6 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

# run the slurm script for blobplot
sbatch blobDB.slurm \
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/purging/aligned.bam \
Kely_original_assembly_blobDB \
1511202 \
/blue/kawahara/yimingweng/Kely_genome_project/blobplot/kely_hifisam_default.megablast.out \
Kely_original_assembly_blobDB

# 6 arguments
# first: the full path to the input fasta of the assembly
# second: the full path to the bam file where the reads were mapped to the target assembly
# thrid: the output prefix for the blob database which will be a directory
# forth: taxid, check NCBI for the taxid for your species
# fifth: the metablast out put from the blastn command (usually the output from previous step)
# sixth: the output prefix for the final table in tsv format which will be in the current directory

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=assembly_blobDB
#SBATCH -o blobDB.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 8:00:00
#SBATCH -c 8

module load blobtools/2.2

fasta=${1}
bam=${2}
dbprefix=${3}
taxid=${4}
metablast=${5}
tsvprefix=${6}

# create blobDB
blobtools create \
--fasta ${fasta} \
--cov ${bam} \
${dbprefix}

# get the node.dmp file from NCBI, and put it in ./Kely_purge_DB
cd ./${dbprefix}
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar xzf new_taxdump.tar.gz
path=$(pwd)
cd ..

# add the tax and blast results
blobtools add \
--taxdump ${path} \
--taxid ${taxid} \
--hits ${metablast} \
${dbprefix}

# create summary tsv file
blobtools filter --table ${tsvprefix}.tsv ${path}
########################################################################
```
- Use R to visualize the result
```
###### R environment ######
# this script is used to plot the blobplot for Keiferia genome
# please modify this script for the use on the other cases
library(ggplot2)

path <- getwd()
contig_dat <- read.table(paste0(path, "/Kely_original_assembly_blobDB.tsv"), header=T, sep="\t")

x11()
ggplot(contig_dat, aes(x=gc, y=log(aligned_cov), size=length, col=bestsumorder_phylum)) + 
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_point(alpha=0.2) +
  scale_size(range = c(2, 10), name="length") +
  scale_x_continuous(limits=c(0, 1)) +
  geom_hline(yintercept=log(120), linetype="dashed", 
             color = "red", size=0.5) +
  geom_hline(yintercept=log(5), linetype="dashed", 
             color = "red", size=0.5) +
  ylab("log(Contig coverage)") +
  xlab("GC proportion") +
  guides(size = "none") +
  guides(color=guide_legend(override.aes = list(size=3), title="best blast phylum")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
```
- The result is showing here. I have put some notes here. The top two contigs with highest coverage are mitochondrial genome (by blasting the contig to get the result). And purging step has removed another 3 contigs that belong to Arthropod due to the extreme high (>120X) and low (<5X) coverage. 
- And because we don't know whether the 3 contigs being trimmed are false or true duplication, I believe **both genome assemblies should be published together** to ensure that all the bioinformation is kept.

<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/blobplot/Kely_origianl_assembly_blobplot_modified.jpg?raw=true?raw=true">

<br />

## **09/07/2022**  
**\# busco**
**\# RepeatModeler2**  

1. Because there is a putative non-target contig in the assembly (ptg000079l), check the BUSCO result to see if any contribution was made from this contig:
- If so, them remove the contig and redo the BUSCO
- if not, them remove the contig and move on to annotation (no need to repeat the BUSCO as the result should not change)

```
[yimingweng@login5 run_endopterygota_odb10]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/busco/Kely_purge_5X_120X/Kely_purge_5X_120X/run_endopterygota_odb10

[yimingweng@login5 run_endopterygota_odb10]$ cat full_table.tsv | grep "ptg000079l"
# not thing has returned, ther is no busco genes related to this contig
```

2. manually make the final assembly by removing the putative non-target contig from the purged genome
```
# make final genome assembly
[yimingweng@login5 kely_final]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final

sed -e '/ptg000079l/,+1d' /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_purging/Kely_purge_5X_120X.fasta > kely_final_assembly.fasta

cat /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta | grep -A1  "ptg000073c" > kely_mito_genome.fasta
```

3. With the final genome ready, let's start to annotate this genome
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

sbatch repeatmodeler.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta kely_repeatmodeler

###########################  script content  ###########################
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
cat kely_genome-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > kely-genome-families.prefix.fa.known
cat kely_genome-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > kely-genome-families.prefix.fa.unknown
```

## **09/07/2022**  
**\# RepeatMasker**  

With the run on repeatmodeler2 finished, let's move on to masking the genome with the information generated from repeatmodeler2. However, I think it is good idea to follow [Dr. Card's suggestions](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/) about how to comprehensively mask a genome with repeat sequences. There are four steps, and I am going to do them in four separated scripts.
1. Step1: mask the simple repeats
```
sbatch kely_repeatmask_step1.slurm

###########################  script content  ###########################
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
########################################################################
```

2. Step2: mask repeats based on existing databases by specifying target taxa (here I use Lepidoptera)
<span style="color:red"> **(this script has been modified we rerun, see below)** <span>

```
sbatch kely_repeatmask_step2.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step2
#SBATCH -o kely_repeatmask_step2.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step2

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-e RMBlast \
-gff \
-no_is \
-species Lepidoptera \
-dir kely_repeatmasker_step2 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step1/kely_final_assembly.fasta.masked &> kely_repeatmasker_step2.out

mv kely_repeatmask_step2.out kely_repeatmasker_step2
########################################################################
```

3. Step3: mask repeats based on the known repeats from RepeatModeler
<span style="color:red"> **(this script has been modified we rerun, see below)** <span>

```
sbatch kely_repeatmask_step3.slurm

###########################  script content  ###########################
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
-lib /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely-genome-families.prefix.fa.known \
-dir kely_repeatmasker_step3 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step2/kely_final_assembly.fasta.masked.masked &> kely_repeatmasker_step3.out

mv kely_repeatmask_step3.out kely_repeatmasker_step3
########################################################################
```

4. Step4: mask repeats based on the unknown repeats from RepeatModeler
<span style="color:red"> **(this script has been modified we rerun, see below)** <span>

```
sbatch kely_repeatmask_step4.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step4
#SBATCH -o kely_repeatmask_step4.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step4

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-e RMBlast \
-gff \
-no_is \
-lib /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely-genome-families.prefix.fa.unknown \
-dir kely_repeatmasker_step4 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step3/kely_final_assembly.fasta.masked.masked.masked &> kely_repeatmasker_step4.out

mv kely_repeatmask_step4.out kely_repeatmasker_step4
########################################################################
```
<br />

## **09/21/2022**  
**\# RepeatMasker** (rerun)

After checking the repeat masking result, I found that all the repeats were masked by hardmask (ATCG -> N), which I think is too strict. Let's check where we loss the softmasked nucleotide.
```
[yimingweng@login5 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation
cat kely_repeatmasker_step1/kely_final_assembly.fasta.masked | grep -v ">" | tr -dc a-z | wc -c
# 6340932

cat kely_repeatmasker_step2/kely_final_assembly.fasta.masked.masked | grep -v ">" | tr -dc a-z | wc -c
0
```

The softmask nucleotides were gone in the second step where I applied the lepidopteran library. So let's rerun RepeatMasker again from step 2 with -xsmall option applied
- rerun Step2: mask repeats based on existing databases by specifying target taxa (here I use Lepidoptera)
```
sbatch kely_repeatmask_step2.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step2
#SBATCH -o kely_repeatmask_step2.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step2

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-species Lepidoptera \
-dir kely_repeatmasker_step2 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step1/kely_final_assembly.fasta.masked &> kely_repeatmasker_step2.out

mv kely_repeatmask_step2.out kely_repeatmasker_step2
########################################################################
```
- track the number of softmasking nucleotide
```
cat kely_repeatmasker_step2/kely_final_assembly.fasta.masked.masked | grep -v ">" | tr -dc a-z | wc -c
6416179
```


- rerun Step3: mask repeats based on the known repeats from RepeatModeler
```
sbatch kely_repeatmask_step3.slurm

###########################  script content  ###########################
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
-xsmall \
-e RMBlast \
-gff \
-no_is \
-lib /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely-genome-families.prefix.fa.known \
-dir kely_repeatmasker_step3 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step2/kely_final_assembly.fasta.masked.masked &> kely_repeatmasker_step3.out

mv kely_repeatmask_step3.out kely_repeatmasker_step3
########################################################################
```

- track the number of softmasking nucleotide
```
cat kely_repeatmasker_step3/kely_final_assembly.fasta.masked.masked.masked | grep -v ">" | tr -dc a-z | wc -c
182220265
```

- rerun Step4: mask repeats based on the unknown repeats from RepeatModeler
```
sbatch kely_repeatmask_step4.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=kely_repeatmask_step4
#SBATCH -o kely_repeatmask_step4.out
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 16

mkdir kely_repeatmasker_step4

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1

RepeatMasker -pa 16 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-lib /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely-genome-families.prefix.fa.unknown \
-dir kely_repeatmasker_step4 \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step3/kely_final_assembly.fasta.masked.masked.masked &> kely_repeatmasker_step4.out

mv kely_repeatmask_step4.out kely_repeatmasker_step4
cp kely_repeatmasker_step4/kely_final_assembly.fasta.masked.masked.masked.masked ./
mv kely_final_assembly.fasta.masked.masked.masked.masked kely_final_masked.fasta
########################################################################
```


- track the number of softmasking nucleotide
```
cat kely_repeatmasker_step3/kely_final_assembly.fasta.masked.masked.masked | grep -v ">" | tr -dc a-z | wc -c

```

- store the final masked genome in the assembly folder (`/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/`)
```
cp kely_repeatmasker_step4/kely_final_assembly.fasta.masked.masked.masked.masked /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/
mv /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_assembly.fasta.masked.masked.masked.masked /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_repeatmasked.fasta
```
<br />

## **09/22/2022**  
**\# genome annotation**  
**\# BRAKER2** 

With the genome with repeats masked, let's try annotate the genome with functional genes using [BRAKER2](https://github.com/Gaius-Augustus/BRAKER#braker-with-proteins-of-any-evolutionary-distance). For this species, I will use the pipeline to run braker with proteins of any evolutionary distance .Later to use the RNA sequence of *Tuta absoluta* (species in the same subfamily to *Keferia*) published by [Camargo et al, 2015](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1841-5) to improve the gene model. <span style="color:red"> (at this point, I am not very certain whether including the RNA sequence data from *Tuta* will improve or disprove the gene model) </span>. I also took suggestions from this [instruction](https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_to_Braker2.html#gsc.tab=0) by Dr. Masonbrink for running bracker2.

1. step1. get the protein database from OrthoDB which is recommended by the authors
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

# download the protein files and put them together
wget https://v101.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
tar xvf odb10_arthropoda_fasta.tar.gz
cat arthropoda/Rawdata/* > proteins.fasta
```

2. step2: Running [ProtHint](https://github.com/gatech-genemark/ProtHint#protein-database-preparation) to get the protein gff file. The gff file here is the tab delimited format file containing the coordinates and related information of the identified proteins based on the provided database (protein.fasta in this case).
```
sbatch -J Kely_ProtHint ProtHint.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_final_masked.fasta /blue/kawahara/yimingweng/Kely_genome_project/annotation/proteins.fasta 

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
dates;hostname;pwd

genome=${1}
protein=${2}

module load prothint/2.6.0

prothint.py --threads ${SLURM_CPUS_ON_NODE:-1} ${genome} ${protein}
########################################################################
```
3. run Braker2 with the protein database

```
[yimingweng@login5 braker2]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein

sbatch -J kely_braker2_protein kely_braker2_protein.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_repeatmasked.fasta /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/prothint_augustus.gff keiferia_lycopersicella

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
protein_gff=${2}
species=${3}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} --hints=${protein_gff} --softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
########################################################################
```

## **10/03/2022**  
**\# genome annotation**  
**\# BRAKER2** 

The Braker annotation with protein database is done. Now use the RNA sequence as reference to rerun Braker and later I will have to merge the two result.

1. download Tuta RNA sequence files (6 files for six different life stages)
```
[yimingweng@login6 RNA_seq_data]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/RNA_seq_data

sbatch SRR_download.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=fastq_dump    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yimingweng@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=04:00:00               # Time limit hrs:min:sec
#SBATCH --output=fastq_dump.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.9
fastq-dump --split-files --gzip SRR2147324 SRR2147322 SRR2147323 SRR2147321 SRR2147320 SRR2147319
########################################################################
```
Note that the sample of SRR2147324 is single end read, not paired end read as the authors claimed in the paper.


2. clean up the reads with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

```
[yimingweng@login5 RNA_seq_data]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/RNA_seq_data

sbatch -J kely_trimmomatic trimmomatic.slurm 

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load trimmomatic/0.39

for sample in $(ls *fastq.gz | cut -d "_" -f 1 | sort | uniq)
do
    fq1=$(ls ${sample}_1*)
    fq2=$(ls ${sample}_2*)
    trimmomatic PE -threads 16 \
    ${fq1} ${fq2} \
    ${sample}_1_clean.fq.gz ${sample}_1_unpaired.fq.gz \
    ${sample}_2_clean.fq.gz ${sample}_2_unpaired.fq.gz LEADING:3 TRAILING:3 MINLEN:36
done
########################################################################
```
Note that the single end sample "SRR2147324" needs to be run by itself because the code was written for paired end reads

## **10/04/2022**  
**\# genome annotation**  
**\# BRAKER2** 

With the clean RNA sequence data, now I am going to map these reads to the genome to create bam files. The bam files will be used as an input to braker2 to build the gene model. However, as I've mentioned above, <span style="color:red">  I am not very certain whether including the RNA sequence data from *Tuta* will improve or disprove the gene model)<span>.

1. Map the reads to the genome assembly  <span style="color:red"> (Note: the fastq files had been deleted after I got the sorted bam files in step3) <span>.
```
[yimingweng@login5 mapped_bams]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams

sbatch -J kely_hisat kely_hisat.slurm 

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load hisat2/2.2.1-3n

hisat2-build /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta.base

for sample in $(ls /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/RNA_seq_data/*fastq.gz | cut -d "_" -f 1,2,3,4,5,6,7 | sort | uniq)
do
    fq1=$(ls ${sample}_1_clean*)
    fq2=$(ls ${sample}_2_clean*)
    name=$(echo ${sample} | cut -d "/" -f 10 | cut -d "_" -f 1)
    hisat2 -p 32 \
    -x /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta.base \
    -1 ${fq1} -2 ${fq2} -S ${name} --phred33 --novel-splicesite-outfile ${name}.junctions --rna-strandness FR
done

hisat2 -p 32 -x /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta.base \
-U /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/RNA_seq_data/single_end/SRR2147324_1_clean.fq.gz \
-S SRR2147324 --phred33 --novel-splicesite-outfile SRR2147324.junctions --rna-strandness F
########################################################################
```
Note that the mapping rate is quite low, but I still have quite a bit of reads being mapped to Keiferia's genome
```
[yimingweng@login5 mapped_bams]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams

for sam in $(ls SRR* | cut -d "." -f 1 | cut -d "_" -f 1 | sort | uniq)
do
    all=$(cat ${sam} | wc -l)
    mapped=$(cat ${sam} | grep "ptg" | wc -l)
    unmapped=$(cat ${sam} | grep \*$'\t'0$'\t0' | wc -l)
    maprate=$(echo ${mapped}/${all} | bc -l)
    unmaprate=$(echo ${unmapped}/${all} | bc -l)
    echo -e "${sam}\t${all}\t${mapped}\t${maprate}\t${unmapped}\t${unmaprate}"
done
```
| readname | read_count | mapped_read | mapped_rate | unmapped_read | unmapped_rate|
| -------- | ---------- | ----------- | ----------- | ------------- | ------------ |
SRR2147319   | 36017531 | 8323571 | 0.231  | 27693958 | 0.769 |
SRR2147320   | 45016275 | 9279237 | 0.206  | 35737036 | 0.794 |
SRR2147321   | 37705161 | 7111755 | 0.189 | 30593404  | 0.811 |
SRR2147322   | 36685228 | 5577554 | 0.152 |  31107672 | 0.848 |
SRR2147323   | 56770013 | 12247733| 0.216 |  44522278 | 0.784 |
SRR2147324   | 19012633 | 2713121 | 0.143 |  16299512 | 0.857 |
|
2. Convert the sam files to bam files (here I forgot to assign the extension for the outputs from hisat2 so the code looks funny). 

```
[yimingweng@login5 mapped_bams]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams

sbatch -J sam2bam sam2bam.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load samtools/1.15

for sam in $(ls SRR* | cut -d "." -f 1 | sort | uniq)
do
    samtools view --threads 16 -b -o ${sam}.bam ${sam}
    samtools sort -m 7G -o ${sam}_sorted.bam -T ${sam}_temp --threads 16 ${sam}.bam
done
########################################################################
```

3. Remove redundant fastq, sam and bam files (just keep the sorted bam files), and run braker with the sorted bam files.
```
[yimingweng@login5 tuta_RNA_seq]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq

sbatch -J kely_braker2_RNA kely_braker2_RNA.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_repeatmasked.fasta keiferia_lycopersicella_RNA

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
species=${2}

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} \
--bam=/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147319_sorted.bam,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147320_sorted.bam,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147321_sorted.bam,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147322_sorted.bam,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147323_sorted.bam,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/mapped_bams/SRR2147324_sorted.bam \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
########################################################################

```
4. To evaluate the completeness of the gene model, I can run busco on the generated amino acid sequences. I will apply this approach to the models from 1) lepidoptera protein model, 2) Tuta RNA-sequence model, 3) combined model,and 4) combined model with Tuta transcriptome refining.
- busco evaluation on gene model based on lepidoptera protein
```
[yimingweng@login5 lep_protein]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein

sbatch lep_protein_model_busco.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_lep_protein_gene_model_busco
#SBATCH -o Kely_lep_protein_gene_model_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/braker/augustus.hints.aa \
 -o ./busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m protein -c 12
########################################################################
# Result: C:95.4%[S:86.3%,D:9.1%],F:2.6%,M:2.0%,n:2124
```
- busco evaluation on gene model based on *Tuta* RNA sequences
```
[yimingweng@login5 tuta_RNA_seq]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq

sbatch tuta_rna_model_busco.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_RNA_seq_gene_model_busco
#SBATCH -o Kely_RNA_seq_gene_model_busco.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 5:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/braker/augustus.hints.aa \
 -o ./busco_out \
 -l /data/reference/busco/v5/lineages/endopterygota_odb10 \
 -m protein -c 12
########################################################################
# Result: C:92.2%[S:86.1%,D:6.1%],F:4.2%,M:3.6%,n:2124
```

## **10/04/2022**  
**\# TSEBRA**  

With the two lines of external evidence for gene prediction, we can integrate those outputs (predictions from protein database and RNAseq reads) to get the final gene prediction model using [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA). 

```
[yimingweng@login1 tsebra]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra

sbatch Kely_TSEBRA.slurm


###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_TSEBRA
#SBATCH -o Kely_TSEBRA.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 24:00:00
#SBATCH -c 24

module load python3

/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/TSEBRA/bin/tsebra.py \
--keep_gtf /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/braker/augustus.hints.gtf,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/braker/augustus.hints.gtf \
-c /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/kely_tsebra.cfg \
-e /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/braker/hintsfile.gff,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/braker/hintsfile.gff \
-o kely_protein_rnaseq_conbine.gtf


/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl \
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_repeatmasked.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/kely_protein_rnaseq_conbine.gtf \
kely_braker_final_aa.fa
########################################################################
```

- busco evaluation on the final merged gene model from TSEBRA
```
[yimingweng@login5 tuta_RNA_seq]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra

sbatch braker_final_model_busco.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=braker_final_model_busco
#SBATCH -o braker_final_model_busco.slurm.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 10:00:00
#SBATCH -c 12

# define configure file for BUSCO and augustus
# For augustus, if encounter an authorization issue (error pops up when running busco), try to download the augustus repo and use its config dir
export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0
module load hmmer/3.2.1

# run busco command
busco -f -i /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/kely_braker_final_aa.fa \
 -o ./busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 12
########################################################################
# Result: C:81.0%[S:75.1%,D:5.9%],F:4.5%,M:14.5%,n:2124
```
The results don't show any improvement based on the BUSCO completeness. I have issued this to TSEBRA [here](https://github.com/Gaius-Augustus/TSEBRA/issues/23) and it seems that there are a lot of transcripts that aren't supported by hints from the extrinsic evidence so they have been filtered out by TSEBRA. So I am going to try two things:
1. take the suggestions from the author (replied in the issue), to run the tsebra with `--keep_gtf` for the protein prediction and `--gtf` for the RNA-seq prediction, plus to tune the parameter `intron_support` to 0.5 to allow more transcripts with low evidence support to pass the filter.

```
[yimingweng@login1 tsebra]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra

sbatch Kely_TSEBRA_2.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Kely_TSEBRA2
#SBATCH -o Kely_TSEBRA2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 24:00:00
#SBATCH -c 24

module load python3

/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/TSEBRA/bin/tsebra.py \
--keep_gtf /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/braker/augustus.hints.gtf \
--gtf /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/braker/augustus.hints.gtf \
-c /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/kely_tsebra2.cfg \
-e /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/lep_protein/braker/hintsfile.gff,/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tuta_RNA_seq/braker/hintsfile.gff \
-o kely_tsebra2.gtf


/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl \
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_repeatmasked.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/tsebra/kely_tsebra2.gtf \
kely_braker_tsebra2_aa.fa
########################################################################
# Result: C:95.7%[S:79.0%,D:16.7%],F:2.4%,M:1.9%,n:2124
```
The busco duplication is still high. This is very likely due to the presence of isoforms. For the following busco evaluation, I should use the updated busco script to run the busco for the gene model. And here are the conclusion and concerns that lead to next step:
- conclusion: for Keiferia genome, the use of Tuta RNAseq read only marginally increase the number of predicted gene, but introduced large number of noise (likely false predicted genes), which can be seen by the boosting duplication rate ( 9% -> 16% )
- concern: both prohint and comnined models (prohint+RNAseq) have predicted too many protein coding genes. The general number of coding genes in Lepidotera is ~15k to 18k. But the prohint and combined model has 36k and 42k genes. With the most strict parameters applied to the tsebra to get the most conserved model, I still get 26k genes with low busco completeness (C:81.0%[S:75.1%,D:5.9%],F:4.5%,M:14.5%). I will need to fix this problem to create a proper gene model.  
<br />


## **10/22/2022** 
There is an issue addressed in the braker github, talking about the problem of having too many genes from braker2. Here is the [link](https://github.com/Gaius-Augustus/BRAKER/issues/319) to it. The author suggested to use their scripts to select the genes only fully, or at least partially supported by the hints. The final model will then be the subset of the total output of braker containing only the higher confident genes. To simplified the steps, I merged all the necessary steps into one slurm scripts called **braker_noRNA_anaotation.slurm**. These steps are:
1. take the repeat-masked genome as input, and run prohints to create protein database called prothint.gff
2. use the prohints database to run braker, of course without the RNAseq, to predict the gene model.
3. use "selectSupportedSubsets.py" to extract the genes with full or partial support from the hints, and write it to a new gtf file called "supported_gene.gtf"
4. Use supported_gene.gtf to generate a new protein fasta file called "braker_supported_gene.aa"
5. To run busco without presence of isoform, further subset the braker_supported_gene.aa to a temporal fasta file (will be removed after the process finished) and this temporal fasta will be used to evaluate the completeness of the genome model with only the longest unique isoforms. That means the duplication from this busco result can not be referred to the presence of isoform.

```

sbatch -J kely_noRNA_model2 braker_noRNA_anaotation2.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_assembly_softmasked.fasta kely_rerun2

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_%j
#SBATCH --output=%x_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
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

# step5: run busco on the subset of the final model (remove smaller isofor)
python3 /blue/kawahara/yimingweng/universal_scripts/longest_aa.py < ${species}_anysupport_aa.fa > ${species}_tmp.fasta

busco -f -i ${species}_tmp.fasta \
 -o ./${species}_gene_model_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32

rm ${species}_tmp.fasta
########################################################################
# BUSCO result: C:92.9%[S:91.5%,D:1.4%],F:0.8%,M:6.3%,n:5286
```

## **10/25/2022**  
**\# functional annotation**  
**\# diamond**  
**\# RefSeq non-redundant protein database**  
Because the transcriptome assembly from *Tuta* is not available, I will have to skip PASA gene model refining work and move forward to do functional annotation. To do the functional annotation, I will use [diamond](https://github.com/bbuchfink/diamond), [InterProScan](https://github.com/ebi-pf-team/interproscan), based on the database [InterPro](https://www.ebi.ac.uk/interpro/).

Let's start with diamond, I will do this work following the [tutorial](https://github.com/bbuchfink/diamond/wiki/1.-Tutorial) on its GitHub Wiki. Because diamond will use blastp to annotate the protein sequence from the output of braker2, I will first create the binary diamond database for blasting, then use this database in dnmd format for blasting. These step are written into a single script called "diamond.slurm". And it takes four arguments:  
1. database for blast in fasta or dnmd format
2. protein sequences of the gene model from braker2 in fasta format
3. the e-value cutoff, usually 1*e-5
4. the output prefix for saving the tsv file  

**Let start blast the protein to the databases:**  
1. The RefSeq non-redundant proteins. This is a huge database containing all the proteins from specific group.
- I use `wget` to get it and it needs to be `gunzip` before use:
```
[yimingweng@login1 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
```
- Run diamond to get the functional annotation from the nr database. <span style="color:red"> Note: this script has been updated, see script written on 10/27/2022. </span>
```
[yimingweng@login1 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond

sbatch -J kely_nr /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond/nr.dmnd /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa 0.00001 kely_nr_k5_1e5

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=kely_diamond
#SBATCH -o kely_diamond.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 48:00:00
#SBATCH -c 24

module load diamond/2.0.9

# example

database=${1} # full path to the database in dmnd format, the diamond will use it to find the function for the querying gene model
gene_model=${2} # gene model for functional annotation in fasta format
cutoff=${3}
outname=${4}
path=$(echo ${database} | rev | cut -d "/" -f 2- | rev)
database_name=$(echo ${database} | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
dmnd=$(ls ${path}/${database_name}.dmnd)

if [ -z "${dmnd}" ]
then
  echo "converting database from fasta to dmnd format"
  diamond makedb --in ${database} -d nr
else
  echo -e "the database has been converted to dmnd format, skip this step and run diamond"
  diamond blastp -k5 -e ${cutoff} -d ${database} -q ${gene_model} -o ${outname}.tsv
fi
########################################################################
```

2. The uniprot arthropod database (Reviewed Swiss-Prot for arthropod)  <span style="color:red"> Note: this script has been updated, see script written on 10/27/2022. </span>
```
[yimingweng@login1 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond

sbatch -J kely_uniprot /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond/uniprot_arthropod.dmnd /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa 0.00001 kely_uniprot_k5_1e5
```

## **10/27/2022**  

The outputs from Diamond are in pretty good shape but the standard annotation format is usually gtf or gff3, so I would like to convert them to gff3 format using the python script called `blast2gff.py` from [genomeGTFtools](https://github.com/wrf/genomeGTFtools). I have put this step into the diamond slurm script.


## **11/02/2022** 
**\# functional annotation**  
**\# InterProScan**  

1. Usually we annotate functions for a gene model we want multiple lines of evidence. Here I am going to use [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/index.html) to annotate the genes

```
[yimingweng@login5 interproscan]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/functional_annotation/interproscan

sbatch -J kely_interproscan /blue/kawahara/yimingweng/universal_scripts/interproscan.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa kely_interproscan

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=%x_interproscan_%j
#SBATCH -o %x_interproscan_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 120:00:00
#SBATCH -c 32

module load iprscan/5.57

protein=${1}
prefix=${2}

cat ${protein} | sed 's/\*//g' > ${prefix}.fasta
interproscan.sh -i ${prefix}.fasta -cpu 32 -f tsv -goterms
########################################################################
```

## **11/04/2022** 
**\# Orthofinder** 

```
[yimingweng@login5 orthofinder]$ pwd
/blue/kawahara/yimingweng/gele_genomes_analyses/orthofinder

sbatch -J gele /blue/kawahara/yimingweng/universal_scripts/orthofinder.slurm /blue/kawahara/yimingweng/gele_genomes_analyses/orthofinder/input_aa

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=%x_orthofinder_%j
#SBATCH -o %x_orthofinder_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 96:00:00
#SBATCH -c 32
#SBATCH --pos=kawahara-b

module load orthofinder/2.5.2

input=${1} # input folder with all amino acid sequences of species we want to include

orthofinder -f ${input}
###########################  script content  ###########################
```










<br />  
<br /> 
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  
<br />  


### Path change/File moving notes
- 10/28/2022: the file `nr.dmnd` was in `/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond` before and now is in `/orange/kawahara/yimingweng/databases/`
- 10/28/2022: the file `uniprot_arthropod.dmnd ` was in `/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond` before and now is in `/orange/kawahara/yimingweng/databases/`