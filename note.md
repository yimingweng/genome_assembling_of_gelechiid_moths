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

3. To further check the read quality, use FASTQC to get a look at the qulality summary.
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
- This is possibly due to the repeat regions at the 100X peak which is odd, as it seems that part of the gemome is **tetraploid**, as the repeats are ~25/~50/~100.
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
3. The result of busco looks good, it has the completeness (C) of **98.1%** but like it is expected, the duplication rate is a bit high (it's **2.3%** and usually I see duplucation rate between 0.5%-2%). So let's go over the [purge_haplotigs pipeline](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/). There are several steps in this pipeline, let's do it step by step.
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
-o Kely_purge_15X_120X
########################################################################
```

- The resulting genome assembly has only **62 contigs** (originally it has 116 contigs), but the assembled size doesn't reduce much (**455 Mbps >> 444 Mbps**). It seems that it has removed most of the smallest contigs even though the cutoffs were not about size but coverages.

4. Rerun BUSCO on the purged genome.
```
[yimingweng@login6 Kely_purge_15X_120X]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/busco/Kely_purge_15X_120X

sbatch universal_run_busco.slurm /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_15X_120X.fasta Kely_purge_15X_120X

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

# copy the contig files from `/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/`
cp -r /blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus ./
cp -r /blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/config.ini ./

export BUSCO_CONFIG_FILE="blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/config.ini"
export AUGUSTUS_CONFIG_PATH="/blue/kawahara/yimingweng/Kely_genome_project/busco/kely_hifisam_default/Augustus/config/"

# load busco, make sure this is the latest version
module load busco/5.3.0

mkdir ${outprefix}

# run busco command
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
- create the database (blastdb) for [hits](https://blobtools.readme.io/docs/taxonomy-file) file, one of the requred input to make blobplot.
```
[yimingweng@login2 blastdb]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot/blastdb

# download the perl script (update_blastdb.pl) from NCBI
# make it excusible
[yimingweng@login2 blastdb]$ chmod +x  update_blastdb.pl

# to make it work, make sure you load the perl module
module load perl/5.24.1

# test this perl script
./update_blastdb.pl --help
```
