# Experiment note for Phthorimaea absoluta genome assembly 
## 
The note is to record the steps we did for assembling a genome for *Phthorimaea absoluta*. We started from a relatively small hifi read data which has only ~2.2 million reads, covering about 9-11X of the genome, assuming the genome size is around 500-600Mbp. Although we have second sequence data coming later, it was even smaller, about 0.5 million reads so it didn't help improving the assembly. Because the main issue when assembling a genome with low coverage hifi long reads is the high duplication from the duplicated haplotigs, purging the assembly without losing the read contigs is the goal. So here this note we started from purging the genome assembly.

## **10/01/2022**
I merged the two read sets, one is 9X and the other is 1-2X. I tried to get the assembly from the pooled reads. Once I got the assembly, I used haplotig purging pipeline to remove the duplicated contigs. Because multiple runs of this pipeline can potentially improve the result, I tried to run this pipeline 3 times.
1. first run 
```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_merge_purge_step1
#SBATCH -o Phthorimaea_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/hifiasm_default/Phthorimaea_merge_default.asm.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/Phthorimaea_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/aligned.bam  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/hifiasm_default/Phthorimaea_merge_default.asm.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge_cutoff
#SBATCH -o Phthorimaea_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o Phthorimaea_coverage_stats.csv \
-j 80 \
-s 80
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge_haplotigs
#SBATCH -o Phthorimaea_purge_haplotigs.log
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
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/hifiasm_default/Phthorimaea_merge_default.asm.fasta  \
-c /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/Phthorimaea_coverage_stats.csv \
-o Phthorimaea_purge_2X_50X
########################################################################
```


2. second run 

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_merge_purge_step1
#SBATCH -o Phthorimaea_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run2/Phthorimaea_purge_step2_2X_20X.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/Phthorimaea_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run3/aligned.bam  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run2/Phthorimaea_purge_step2_2X_20X.fasta
########################################################################
```


```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge_cutoff_run2
#SBATCH -o Phthorimaea_purge_cutoff_run2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run3/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o Phthorimaea_coverage_stats_run2.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run2/Phthorimaea_purge_step2_2X_20X.fasta  \
-c /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run3/Phthorimaea_coverage_stats_run3.csv \
-o Phthorimaea_purge_step2_2X_20X
########################################################################
```
3. third run

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_merge_purge_step1
#SBATCH -o Phthorimaea_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/Phthorimaea_purge_2X_20X.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/Phthorimaea_merged.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run2/aligned.bam  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/Phthorimaea_purge_2X_20X.fasta
########################################################################
```


```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge_cutoff_run3
#SBATCH -o Phthorimaea_purge_cutoff_run3.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run3/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o Phthorimaea_coverage_stats_run3.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run2/Phthorimaea_purge_step2_2X_20X.fasta  \
-c /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge/purge/run3/Phthorimaea_coverage_stats_run3.csv \
-o Phthorimaea_purge_step2_2X_20X
########################################################################
```
The assembly after 3 runs of the purging pipeline still contains high duplication rate, according to the busco single copy gene evaluation.



## **10/10/2022**
Now we should give up the pooled reads, and run purge pipeline for the old single sample assembly (the one with 9X coverage).
```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge_step1
#SBATCH -o Phthorimaea_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load minimap/2.21
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

minimap2 -t 4 -ax map-pb /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/assembly1/Phthorimaea_hifisam_reassemble_ctg.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/assembly1/m64219e_220604_102704.hifi_reads.fastq.gz \
--secondary=no \
| samtools sort -m 1G -o aligned.bam -T tmp.ali

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/purge/aligned.bam  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/assembly1/Phthorimaea_hifisam_reassemble_ctg.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_original_purge_cutoff
#SBATCH -o Phthorimaea_original_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/purge/aligned.bam.gencov  \
-l 2  \
-m 33  \
-h 20  \
-o Phthorimaea_original_coverage_stats_run3.csv \
-j 80 \
-s 80

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/assembly1/Phthorimaea_hifisam_reassemble_ctg.fasta  \
-c /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble/purge/Phthorimaea_original_coverage_stats_run3.csv \
-o Phthorimaea_purge_step2_2X_20X
########################################################################
```


## **10/19/2022**  
Since it's been long time without updating the progress here, I'll just quickly summarize the recent findings here.
1. The Phthorimaea sequences from the 2 samples are neither good. The first sample has ~9-10X read coverage and the second one has only 1-2X coverage.
2. I've tried different combinations of assembling methods but non of them works perfectly. I will have to lower the standard a little bit and see what the best result I can get.
<img src="https://drive.google.com/uc?export=view&id=19HOXfstmXAP083iLs_k8vkzdXNH9PeIG">

## **10/20/2022**
**\# bwa mem**  
**\# purging**  
Because there was no perfect results from the pooled read assemblies, I would like to go back to the 9X read coverage data. At least I don't need to worry about the heterozygosity issue. Shashank used to assembled the first version of the genome with the 9X dataset and the purging setting was set to 2 (`-l 2`), and it resulted in 97% overall BUSCO but with 79% duplication rate. Then I tried rerunning hifiasm with `-l 3` to increase the power to purge the haplotigs but unfortunately I got overall BUSCO only 84%, even though the duplication rate decreased to 11%, and 3% after further purging.  
One thing I would like to try is to mapped the published short read data by [Tabuloc et al, 2019](https://link.springer.com/article/10.1007/s10340-019-01116-6) which has ~72X short read coverage, to the first version genome (the one with 97% BUSCO but high duplication rate). The reason I'd like to try this is because original step we purge the haplotigs is to map the reads we used to assemble the genome and see the coverage distribution and similarity of the contigs to determine the false splitting of haplotigs by the heterozygosity. However, since we only have 9X read depth so it is likely we don't have power to distinguish the false contigs from the alternative haplotype contigs. So let's try mapping the deep sequence reads to the genome and see if we can make the coverage cutoffs more clear. 
1. download the read files
```
[yimingweng@login2 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads

sbatch  fastq_dump.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=test_fastq_dump    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yimingweng@ufl.edu     # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=fastq_dump_test.log   # Standard output and error log
pwd; hostname; date

module load sra/2.10.9
fastq-dump --split-files --gzip  SRR8676205
########################################################################
```

2. run bwa to map the reads to the genome
```
[yimingweng@login2 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads

sbatch Phthorimaea_bwa.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=Phthorimaea_bwa
#SBATCH --output=Phthorimaea_bwa.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=180:00:00
#SBATCH -c 64

module load bwa/0.7.17
module loead samtools/1.15

bwa index /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/T_absoluta_hifiasm_06_14_2022.asm.bp.p_ctg.fa

bwa mem -t 64 /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/T_absoluta_hifiasm_06_14_2022.asm.bp.p_ctg.fa /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/SRR8676205_2.fastq.gz /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/SRR8676205_3.fastq.gz > Phthorimaea_short_read_aln.sam

samtools view -S -b Phthorimaea_short_read_aln.sam > Phthorimaea_short_read_aln.bam
samtools sort -m 1G -o Phthorimaea_short_read_aln.bam -T tmp.ali
########################################################################
```

## **10/31/2022** 
**\# purging** 

Once we have the mapped bam file, use it to run the haplotig purging again
```
[yimingweng@login5 run2_with_short_reads]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads

sbatch Phthorimaea_purge2_step2.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_purge2_step2
#SBATCH -o Phthorimaea_purge2_step2.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs  hist  \
-b /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_short_read_aln_sorted.bam \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/Phthorimaea_purge.fasta
########################################################################
```

```
###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_original_purge_cutoff
#SBATCH -o Phthorimaea_original_purge_cutoff.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 1:00:00
#SBATCH -c 4

module load bedtools/2.30.0
module load samtools/1.15
module load purge_haplotigs/1.1.2
module load libssl/1.0.2l

purge_haplotigs cov \
-i /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_short_read_aln_sorted.bam.gencov  \
-l 30  \
-m 85  \
-h 190  \
-o Phthorimaea_run2_stats.csv \
-j 95 \
-s 95

purge_haplotigs purge  \
-g /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/Phthorimaea_purge.fasta  \
-c /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_run2_stats.csv \
-o Phthorimaea_purge_run2
########################################################################
Result: C:96.2%[S:82.5%,D:13.7%],F:0.5%,M:3.3%,n:5286
```


## **11/02/2022** 
**\# best assembly model**  
**\# blobplot**  

By considering the BUSCO score, the best assembly I can get is to use the 9X genomic data, and run hifiasm with purging aggressiveness to be `-l 2`, and let the haplotig-purging pipeline to remove the duplications. So here is the final assembly for annotation:`/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta`

1. Check the assembly statistics:
```
[yimingweng@login5 purge]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge

module load python3

python /blue/kawahara/yimingweng/universal_scripts/assemblystats.py /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta
{
  "Contig Stats": {
    "L10": 12,
    "L20": 29,
    "L30": 51,
    "L40": 79,
    "L50": 115,
    "N10": 4435971,
    "N20": 3284839,
    "N30": 2611967,
    "N40": 2067342,
    "N50": 1614219,
    "gc_content": 38.454626166301814,
    "longest": 7993428,
    "mean": 948696.4491279069,
    "median": 630478.0,
    "sequence_count": 688,
    "shortest": 3255,
    "total_bps": 652703157
```
-> L50: **1.61M**  
-> Assembled size" **65 Mbps**

2. run blobplot to identify the potential foreign sequences
    - run megablast to assign the best hit gene to the contig
```
[yimingweng@login5 blobplot]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/blobplot

sbatch -J Phthorimaea_genome /blue/kawahara/yimingweng/universal_scripts/megablast_nt.slurm /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta Phthorimaea_genome
```

3. Store this "final version" of the assembly in `/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies` and name it as "**Phthorimaea_final_assembly.fasta**".
```
yimingweng@login2 assemblies]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies

cp /blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_sho
rt_reads/Phthorimaea_purge_run2.fasta  ./

mv Phthorimaea_purge_run2.fasta Phthorimaea_final_assembly.fasta
```




```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
. blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp


sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/minimap.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta \
/blue/kawahara/shashankp/Phthorimaeaabsoluta_NS2707/Phthorimaea_fastqc/m64219e_220604_102704.hifi_reads.fasta.gz \
Phthorimaea_blobplot


sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/blobplot.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/blobplot/Phthorimaea_blobplot.bam \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/blobplot/Phthorimaea_genome.nt.mts1.hsp1.1e25.megablast.out \
Phthorimaea_blobplot

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/repeatmodeler2.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea_purge_run2.fasta \
Phthorimaea_repeat
```
<img src="https://github.com/yimingweng/Phthorimaea_genome_project/blob/main/blobplot_results/Phthorimaea_blobplot.png">


## **11/07/2022** 
**\# annotation**  
**\# repeatmodeler**  
**\# repeatmasker**   

1. run Repeatmodeler
```
[yimingweng@login2 repeatmodeler]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/repeatmodeler

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/repeatmodeler2.slurm /blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta Phthorimaea_repeat
```

2. run repeatmasker
```
[yimingweng@login2 repeatmasker]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/repeatmasker

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/repeatmakser.slurm /blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/repeatmodeler/Phthorimaea_repeat-families.fa \
Phthorimaea_repeatmasker

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_repeatmask_%j
#SBATCH -o %x_repeatmask_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 120:00:00
#SBATCH -c 32

genome=${1}
modeler_out={2}
prefix=${3}

mkdir ${prefix}

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

module load repeatmasker/4.1.1


name=$(echo ${genome} | rev | cut -d "/" -f 1 | rev | cut -d "." -f 1)
path=$(pwd)

# step 1: mask the simple repeat
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-noint \
-no_is \
-dir ${prefix} \
${genome} &> ./${prefix}/${prefix}_step1.out


# step 2: mask repeats based on existing databases
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-species Lepidoptera \
-dir ${prefix} \
${path}/${prefix}/${name}.masked &> ./${prefix}/${prefix}_step2.out


# step 3: mask genome based on the output of repeatmodeler
RepeatMasker -pa 32 -a -s \
-xsmall \
-e RMBlast \
-gff \
-no_is \
-lib ${modeler_out} \
-dir ${prefix} \
${path}/${prefix}/${name}.masked.masked &> ./${prefix}/${prefix}_step3.out
```
3. check the masing rate
```
[yimingweng@login6 Phthorimaea_repeatmasker]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/repeatmasker/Phthorimaea_repeatmasker

bash /blue/kawahara/yimingweng/universal_scripts/maskrate.sh Phthorimaea_final_assembly.fasta.masked.masked.masked
softmasking rate is 54.40%
hardmasking rate is 0%
```

## **11/07/2022** 
**\# annotation**  
**\# braker2**  
**\# TSEBRA**  

1. run braker2 with the protein database from orthoDB. The prothint data will be used to train the gene model 
```
[yimingweng@login6 prothint]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/prothint

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/braker_prothint.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/repeatmasker/Phthorimaea_repeatmasker/Phthorimaea_final_assembly.fasta.masked.masked.masked \
Phthorimaea_prothint
```

2. run braker2 with the RNAseq data. The reads were download from NCBI originally published by [Camargo et al, 2015](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1841-5). There are 6 different life stages, and one got single end reads and the rest got paired end reads. 
- download the reads
```
[yimingweng@login2 raw_reads]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/raw_reads

sbatch SRR_download.slurm
```
- clean up the reads using [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```
sbatch -J Phthorimaea_trimmomatic trimmomatic.slurm /blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/raw_reads/ /blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/trimmed_reads/

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

inputdir=${1}
outdir=${2}

ls ${inputdir} | grep "fastq.gz" | cut -d "_" -f 1 | sort | uniq >> list.tmp
IFS=$'\n'
for sp in $(cat list.tmp)
do
  count=$(ls ${inputdir} | grep "${sp}" | wc -l)
  echo -e "${sp}\t${count}" >> list2.tmp
done

cat list2.tmp | grep -Pw "2" | cut -d $'\t' -f 1 > pairlist.tmp
cat list2.tmp | grep -Pw "1" | cut -d $'\t' -f 1 > singlelist.tmp


for sample in $(cat pairlist.tmp)
do
    fq1=$(ls ${inputdir}/${sample}_1*)
    fq2=$(ls ${inputdir}/${sample}_2*)
    trimmomatic PE -threads 16 \
    ${fq1} ${fq2} \
    ${outdir}/${sample}_1_clean.fq.gz ${outdir}/${sample}_1_unpaired.fq.gz \
    ${outdir}/${sample}_2_clean.fq.gz ${outdir}/${sample}_2_unpaired.fq.gz LEADING:3 TRAILING:3 MINLEN:36
done

for sample in $(cat singlelist.tmp)
do
    fq1=$(ls ${inputdir}/${sample}*)
    trimmomatic SE -threads 16 \
    ${fq1} \
    ${outdir}/${sample}_1_clean.fq.gz LEADING:3 TRAILING:3 MINLEN:36
done
rm *tmp
########################################################################
```
- map the reads to the genome assembly using [hisat](http://daehwankimlab.github.io/hisat2/manual/)

```
[yimingweng@login2 bam_files]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bam_files

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/hisat.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
Phthorimaea_absoulata \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/trimmed_reads

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_hisat_%j
#SBATCH --output=%x_hisat_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load hisat2/2.2.1-3n

genome=${1}
species=${2}
readdir=${3}

hisat2-build ${genome} ${species}.fasta.base

ls ${readdir} | grep "clean.fq.gz" | cut -d "_" -f 1 | sort | uniq >> list.tmp
IFS=$'\n'
for sp in $(cat list.tmp)
do
  count=$(ls ${readdir} | grep "clean" | grep "${sp}" | wc -l)
  echo -e "${sp}\t${count}" >> list2.tmp
done

cat list2.tmp | grep -Pw "2" | cut -d $'\t' -f 1 > pairlist.tmp
cat list2.tmp | grep -Pw "1" | cut -d $'\t' -f 1 > singlelist.tmp

for sample in $(cat pairlist.tmp)
do
    fq1=$(ls ${readdir}/${sample}_1_clean*)
    fq2=$(ls ${readdir}/${sample}_2_clean*)
    name=$(echo ${sample} | cut -d "_" -f 1)
    hisat2 -p 32 \
    -x ${species}.fasta.base \
    -1 ${fq1} -2 ${fq2} -S ${name}.sam --phred33 --novel-splicesite-outfile ${name}.junctions --rna-strandness FR
done

for sample in $(cat singlelist.tmp)
do
    fq1=$(ls ${readdir}/${sample}_1_clean*)
    name=$(echo ${sample} | cut -d "_" -f 1)
    hisat2 -p 32 \
    -x ${species}.fasta.base \
    -U ${fq1} -S ${name}.sam --phred33 --novel-splicesite-outfile ${name}.junctions --rna-strandness F
done
rm *tmp
rm *junctions
```
- check the mapping rate for those reads
```
[yimingweng@login2 bam_files]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bam_files

for sam in $(ls SRR* | cut -d "_" -f 1 | sort | uniq)
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
SRR2147319   | 40353664 | 36264378 | 0.899  | 4089284  | 0.10  |
SRR2147320   | 50787264 | 46242154 | 0.910  | 4545108  | 0.089 |
SRR2147321   | 43493648 | 40087596 | 0.921  | 3406050  | 0.078 |
SRR2147322   | 41577754 | 37805648 | 0.909  | 3772104  | 0.090 |
SRR2147323   | 64456121 | 58072539 | 0.901  | 6383580  | 0.099 |
SRR2147324   | 22261676 | 17878864 | 0.803  | 4382812  | 0.197 |
|

- convert the sam files to bam
```
[yimingweng@login2 bam_files]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bam_files

sbatch -J sam2bam sam2bam.slurm

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_sam2bam_%j
#SBATCH --output=%x_sam2bam_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load samtools/1.15

for sam in $(ls SRR* | sort | uniq)
do
    name=$(echo ${sam} | cut -d "." -f 1)
    samtools view --threads 16 -b -o ${name}.bam ${sam}
    samtools sort -m 7G -o ${name}_sorted.bam -T ${name}_temp --threads 16 ${name}.bam
done
########################################################################
```

- check the RNA read coverage, the coverage should be sufficient to cover the intros so that the braker can use it to train GeneMark-ET and to predict the gene model with the spliced alignment information in Augustus.
```
[yimingweng@login5 bam_files]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bam_files

sbatch depth.slurm

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=Phthorimaea_merge_purge_step1
#SBATCH -o Phthorimaea_merge_purge_step1.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=8gb
#SBATCH -t 5:00:00
#SBATCH -c 4

module load samtools/1.15

for bam in $(ls *bam)
do
    samtools depth ${bam} | awk '{sum+=$3} END { print "Average = ",sum/NR}' >> depth
done
########################################################################
[yimingweng@login5 bam_files]$ cat depth
Average =  5.5942 (SRR2147319)
Average =  7.25236 (SRR2147320)
Average =  7.34867 (SRR2147321)
Average =  5.9289 (SRR2147322)
Average =  8.33535 (SRR2147323)
Average =  3.22731 (SRR2147324)
```
- remove redundant fastq, sam and bam files (just keep the sorted bam files), and run braker with the sorted bam files
```
[yimingweng@login2 RNAseq]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/braker2_RNA.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
Phthorimaea_absoluta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bamlist

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_braker2_RNA_%j
#SBATCH --output=%x_braker2_RNA_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
dates;hostname;pwd

genome=${1}
species=${2}
bamlist=${3}
bam=$(cat ${bamlist} | sed -z 's/\n/,/g')

module load conda
module load braker/2.1.6

braker.pl \
--AUGUSTUS_CONFIG_PATH=/blue/kawahara/yimingweng/LepidoPhylo_Project/busco_out/Augustus/config \
--genome=${genome} --species ${species} \
--bam=${bam} \
--softmasking --gff3 --cores 32 --AUGUSTUS_ab_initio
########################################################################
```
- Use [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA) to combine the gene models from prothint and RNAseq models.
```
[yimingweng@login2 tsebra]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra

sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/tsebra.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/prothint/braker/braker.gtf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/braker/braker.gtf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.cgf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/prothint/braker/hintsfile.gff \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/braker/hintsfile.gff \
Phthorimaea_tsebra_default

###########################  script content  ###########################
#!/bin/bash

#SBATCH --job-name=%x_tsebra_%j
#SBATCH -o %x_tsebra_%j.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mem-per-cpu=4gb
#SBATCH -t 24:00:00
#SBATCH -c 24

module load python3

genome=${1} # the genome assembly
prothint_gff=${2} # the gtf output from braker using prothint evidence
rnaseq_gff=${3} # the gtf output from braker using RNAseq evidence
config=${4} # the config file, the default cfg can be download from the TSEBRA repo
prothint_hints=${5} # the hintsfile.gff from braker using prothint evidence
rnaseq_hints=${6} # the hintsfile.gff from braker using RNAseq evidence
outprefix=${7} # the name for the output file

/blue/kawahara/yimingweng/universal_scripts/TSEBRA/bin/tsebra.py \
--gtf ${prothint_gff},${rnaseq_gff} \
-c ${config} \
-e ${prothint_hints},${rnaseq_hints} \
-o ${outprefix}.gtf

/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl \
${genome} \
${outprefix}.gtf \
${outprefix}.aa
########################################################################
```

## **11/08/2022**   
**\# gene model checking and trimming**  
- check the number of transcript and busco results in the combined model
```
[yimingweng@login2 tsebra]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra

# check the number of transcript
cat Phthorimaea_tsebra_default.aa | grep ">" | wc -l
# 48412

# run busco on the combined gene model
sbatch -J Phthorimaea /blue/kawahara/yimingweng/universal_scripts/busco_gene_model.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.aa

### Results: C:93.9%[S:75.8%,D:18.1%],F:1.5%,M:4.6%,n:5286
```

- remove the transcripts without hints support, and redeem the non-supported genes with blast hits
```
[yimingweng@login2 tsebra]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra

# trim the transcript with no hint support
sbatch -J Phthorimaea findsupport.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.gtf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/prothint/braker/hintsfile.gff \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
Phthorimaea_absoluta

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_findsupport_%j
#SBATCH --output=%x_findsupport_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

combined_gtf=${1} # the gtf file from the combined model
hintsfile=${2} # the hintsfile.gff from the prothint model
genome=${3}
species=${4}

module load prothint/2.6.0
module load braker/2.1.6
module load busco/5.3.0
module load hmmer/3.2.1
module load diamond/2.0.9

# step 1: extract the genes with full or partial support from the hints
python3 /blue/kawahara/yimingweng/universal_scripts/selectSupportedSubsets.py ${combined_gtf} ${hintsfile} --fullSupport fullsupport.gff --anySupport ${species}_anySupport.gff --noSupport nosupport.gff

# step 2: redeem genes from nosupport.gff, by blasting them to the databases nr and UniProt
# step 2.1: get the amino acid sequence
/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl ${genome} nosupport.gff ${species}_nosupport.aa

# step 2.2: run diamond to blast the nosupport aa sequences to the nr and UniProt 
diamond blastp -k 1 -e 0.00001 -d /orange/kawahara/yimingweng/databases/nr.dmnd -q ${species}_nosupport.aa -f 6 -o nosupport_nr.tsv
diamond blastp -k 1 -e 0.00001 -d /orange/kawahara/yimingweng/databases/uniprot_arthropod.dmnd -q ${species}_nosupport.aa -f 6 -o nosupport_uniprot.tsv
cat nosupport_nr.tsv nosupport_uniprot.tsv | cut -d $'\t' -f 1 | sort | uniq >> transcript_list
IFS=$'\n'
for transcript in $(cat transcript_list); do cat nosupport.gff | grep -Pw ${transcript} | cut -d $'\t' -f 9 | cut -d " " -f 4 | grep -o  "[A-Za-z]_[0-9]*">> redeem_list; done

# step 2.3: add genes with blast hits to the gene model
cat ${species}_anySupport.gff | cut -d $'\t' -f 9 | cut -d " " -f 4 | grep -o  "[A-Za-z]_[0-9]*" >> redeem_list
cat redeem_list | sort | uniq >> redeem_list_uniq
sort -t "_" -k2,2 -n redeem_list_uniq >> redeem_list_uniq_sort
rm redeem_list redeem_list_uniq
IFS=$'\n'
for gene in $(cat redeem_list_uniq_sort); do cat ${combined_gtf} | grep -Pw ${gene} >> ${species}_braker_final.gtf; done

# step 2.3: supported_gene.gtf to generate final gene model fasta to anysupport.aa
/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl ${genome} ${species}_braker_final.gtf ${species}_braker_final.aa

# step3: run busco on the subset of the final model (remove smaller isofor)
python3 /blue/kawahara/yimingweng/universal_scripts/longest_aa.py < ${species}_braker_final.aa > ${species}_tmp.fasta

busco -f -i ${species}_tmp.fasta \
 -o ./${species}_gene_model_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32

rm ${species}_tmp.fasta ${species}_nosupport.aa
########################################################################
# BUSCO result: C:93.1%[S:76.8%,D:16.3%],F:1.5%,M:5.4%,n:5286
```

## **11/17/2022** 
**\# annotation**  
**\# braker2**  
**\# TSEBRA** 

1. Although the model generated from this pipeline looks good in BUSCO (consider the original assembly has 13% duplication rate), the gene number is still way too high. I found that this inflating model is common in braker2 with RNA sequence training model. After trying many different filtering processes including [gFACs](https://gfacs.readthedocs.io/en/latest/), the best model based on gene number and BUSCO is from this pipeline:

<img src="https://github.com/yimingweng/Phthorimaea_genome_project/blob/main/annotation/Phthorimaea_gene_model_pipeline.jpg">

And here this the script combining the filtering and transcribe selection (TSEBRA) steps:
```
[yimingweng@login6 braker2]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra2

sbatch -J Phthorimaea findsupport_tsebra.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/prothint \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra2/Phthorimaea_tsebra_default.cgf

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_findsupport_tsebra_%j
#SBATCH --output=%x_findsupport_tsebra_%j.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

module load busco/5.3.0

genome=${1} # soft masked genome
RNAdeq=${2} # path to the folder containing braker model from RNA seq data
prothint=${3}  # path to the folder containing braker model from protein database
config=${4} # the config file with parameters in there (Here I leave it default)

name=$(echo ${config} | cut -d "." -f 1)

# get hint support model from the RNAseq trained model
python3 /blue/kawahara/yimingweng/universal_scripts/selectSupportedSubsets.py ${RNAdeq}/braker/augustus.hints.gtf ${RNAdeq}/braker/hintsfile.gff --fullSupport rna_fullsupport.gff --anySupport rna_anySupport.gff --noSupport rna_nosupport.gff

# get hint support model from the prothint trained model
python3 /blue/kawahara/yimingweng/universal_scripts/selectSupportedSubsets.py ${prothint}/braker/augustus.hints.gtf ${prothint}/braker/hintsfile.gff --fullSupport pro_fullsupport.gff --anySupport pro_anySupport.gff --noSupport pro_nosupport.gff


/blue/kawahara/yimingweng/universal_scripts/TSEBRA/bin/tsebra.py \
--gtf pro_anySupport.gff,rna_anySupport.gff \
-c ${config} \
-e ${prothint}/braker/hintsfile.gff,${RNAdeq}/braker/hintsfile.gff \
-o ${name}.gtf

/blue/kawahara/yimingweng/universal_scripts/Augustus/scripts/gtf2aa.pl \
${genome} \
${name}.gtf \
${name}.aa

# run busco on the subset of the final model (remove smaller isofor)
python3 /blue/kawahara/yimingweng/universal_scripts/longest_aa.py < ${name}.aa > tmp.fasta

busco -f -i tmp.fasta \
 -o ./single_isoform_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32
rm $tmp.fasta

busco -f -i ${name}.aa \
 -o ./all_isoform_busco_out \
 -l /data/reference/busco/v5/lineages/lepidoptera_odb10 \
 -m protein -c 32
########################################################################
BUSCO for single isoform: C:93.1%[S:76.7%,D:16.4%],F:1.3%,M:5.6%,n:5286
BUSCO for all isoforms: C:93.2%[S:75.9%,D:17.3%],F:1.3%,M:5.5%,n:5286
```

2. With the final model settled down. I would like to look closer to the model and some of the statistics can be used to compared with models from other gelechiid species (Keferia and Phthorimaea, for example). Here I employed [gFACs](https://gfacs.readthedocs.io/en/latest/) again, but not foe the model trimming but for accessing the statistics of the model.

```
[yimingweng@login6 gfacs]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/gfacs

sbatch -J Phthorimaea_final_model /blue/kawahara/yimingweng/universal_scripts/gfacs_statistics.slurm \
braker_2.1.2_gtf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.gtf \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/assemblies/Phthorimaea_final_assembly.fasta \
Phthorimaea_tesbra_gfacs \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/gfacs

###########################  script content  ###########################
#!/bin/bash
#SBATCH --job-name=%x_gfacs_default
#SBATCH --output=%x_gfacs_default.log
#SBATCH --mail-user=yimingweng@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module load perl/5.24.1

format=${1}
model=${2}
genome=${3}
prefix=${4}
outdir=${5}

perl /blue/kawahara/yimingweng/universal_scripts/gFACs/gFACs.pl \
-f ${format} \
-p ${prefix} \
    --statistics \
    --statistics-at-every-step \
--distributions \
    exon_lengths 10 \
    intron_lengths 15 \
    CDS_lengths 20 \
    gene_lengths 100 \
    exon_position \
    exon_position_data \
    intron_position \
    intron_position_data \
--fasta ${genome} \
    --nt-content \
--create-gff3 \
--get-protein-fasta \
--create-gtf \
-O ${outdir} \
${model}
########################################################################
```


## **11/17/2022** 
**\# InterProScan** 

```
[yimingweng@login5 interporscan]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/interporscan

sbatch -J Phthorimaea_interproscan /blue/kawahara/yimingweng/universal_scripts/interproscan.slurm \
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.aa \
Phthorimaea_interproscan
```

## **11/27/2022** 
**\# diamond** 
```
[yimingweng@login6 diamond]$ pwd
/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/diamond

sbatch -J Phthorimaea_nr /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /orange/kawahara/yimingweng/databases/nr.dmnd /blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.aa 0.00001 Phthorimaea_nr_k5_1e5

sbatch -J Phthorimaea_uniprot /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /orange/kawahara/yimingweng/databases/uniprot_arthropod.dmnd /blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/tsebra/Phthorimaea_tsebra_default.aa 0.00001 Phthorimaea_uniprot_k5_1e5
```

## 

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


### **Path change/File moving notes**
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_merge`
- 11/02/2022: rename directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_new_seq/` to be `~/Phthorimaea_second_sequence/`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_second_sequence/genome_size`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_old_hifisam_reassemble`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/*fastq.gz`
- 11/02/2022: remove directory `/blue/kawahara/yimingweng/Phthorimaea_genome_project/Phthorimaea_first_assembly/purge/run2_with_short_reads/Phthorimaea*.bam`
- 11/08/2022: remove files `/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/raw_reads/*.fastq.gz`
- 11/08/2022: remove files `/blue/kawahara/yimingweng/Phthorimaea_genome_project/annotation/braker2/RNAseq/bam_files/*sam*`