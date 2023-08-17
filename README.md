# Kely genome project
The genome assembly project for the Tomato Pinworm, *Keiferia lycopersicella* (Walsingham), *Phthorimaea absoluta* (*Tuta absoluta*, tomato leafminer), and *Scrobipalpa atriplicella* (goosefoot groundling moth).

## Background
- The tomato pinworm is a small moth feeds on mostly tomatoes and sometimes eggplant. The damage is mainly on the leaf through the leaf-mining feeding habit.
- The tomato pinworm is hard to be distinguished from several other closely related species by the external morphology.
- The host range of tomato pinworm is mainly on Solanaceae especially the tomato, but the attacks to other species like egg plant and some weed species in the family are documented.
- A comprehensive document about this species can be found in the [Featured Creatures website](https://entnemdept.ufl.edu/creatures/veg/tomato/tomato_pinworm.htm) published by University of Florida.
- It would be interesting to have a genome assembled, annotated, and published for this species, as it could provide a promising scientific value to understand the pest ecology/evolution and potentially the genes involved in the host plant selection.
- To compare the genome with other gelechiid moth, we also sequenced and assembled the other two species: *Phthorimaea absoluta* and *Scrobipalpa atriplicella*.

<br />

## General Workflows
This section describes the steps we used for a general genome assembling pipeline from raw reads to a well-annotated genomes. The description is mainly for *Keiferia lycopersicella*, but also work for the other two species. Except for the purging step of *Phthorimaea absoluta*. Details of those small modifications please see the note files.
<br />

### **Read Quality Assessment**
The journey starts from a fastq file, which contains the hifi ccs reads from PacBio sequencing platform. First thing to do is the check the quality of these reads, so I used [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to assess the read quality. The script called `hifi_fastqc.slurm` can be found in the scripts folder.
```
[yimingweng@login2 raw_reads]$ sbatch -J kely_fastqc hifi_fastqc.slurm Keiferia_lycopersicella_ccs.fastq.gz
```
The result of the fastqc can be viewed [here](https://htmlpreview.github.io/?https://github.com/yimingweng/Kely_genome_project/blob/main/Keiferia_lycopersicella_ccs_fastqc.html). The GC content peak is a bit sharp but it shouldn't be the problem because the adaptors have been trimmed before sending to us.  
<br />  


### **Genome size and sequence coverage estimation**
The estimations of genome size and sequence coverage (read depth) usually come together because they are dependent to each other. Two steps to get the estimates of these two parameters: kmer-counting and kmer-distribution modeling. The kmer-counting job here is done using [KMC](https://github.com/refresh-bio/KMC), while the kmer-distribution modeling is done by [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/).
1. kmer-count: take the fastq file as input to run KMC. The script called `fastq_KMC.slurm` can be found in the scripts folder. This script takes 3 arguments: (1) input read fastq file with full path; (2) prefix of output file; (3) kmer size, here I use k=21.
```
[yimingweng@login2 genome_size]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/genome_size

sbatch -J kely fastq_KMC.slurm /blue/kawahara/yimingweng/Kely_genome_project/raw_reads/Keiferia_lycopersicella_ccs.fastq.gz Kely 21
```
2. Use [GenomeScope2](http://qb.cshl.edu/genomescope/genomescope2.0/) (online tool) to visualize the kmer distribution and get the estimates of genome size and sequence coverage. The parameters were kept default with kmer size to be 21. The result can be viewed [here](http://qb.cshl.edu/genomescope/genomescope2.0/analysis.php?code=mYfpkD8wYxxcIy48XpO0); and the distribution is shown below:

<img src="http://qb.cshl.edu/genomescope/genomescope2.0/user_data/mYfpkD8wYxxcIy48XpO0/linear_plot.png?raw=true">


- The estimated genome size is about 302 Mbp, since there are peak at 100X, we know there are repeats in the genome and the genome size 302Mbp is probably underestimated.
- The heterozygosity is fairly low, about 1%.
- The mean coverage is about 24.7X, slightly lower than expected but should be enough for de novo assembling.
- I also ran K=27 and got similar result. Details please see the [research note](https://github.com/yimingweng/Kely_genome_project/blob/main/note.md).
<br />  

### **Genome Assembling**
1. The genome assembly was perform using [hifiasm](https://github.com/chhylp123/hifiasm) with default setting. The script called `kely_hifisam_default.slurm` can be found in the scripts folder.
```
[yimingweng@login2 kely_hifisam_default]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/
sbatch kely_hifisam_default.slurm
```

2. Check the assembly statistics using the customized python script called `assemblystats.py` which is also in the scripts folder.
- It takes the assembled genome as first input, it needs to be converted to fasta format. The conversion was done in the `kely_hifisam_default.slurm` script.
- The second argument is the output prefix.
```
[yimingweng@login2 kely_hifisam_default]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default
./assemblystats.py Kely_hifisam.asm.bp.p_ctg.gfa Kely_hifisam_default_stats
```
- The assembly has **116** contigs with N50=**14.12Mbp** and L50=**12**. The assembled size is **455.7Mbp**, slightly larger than the estimate of GenomeScope (318 Mbp).
3. Use [BUSCO](https://busco.ezlab.org/) to assess the completeness of the assembly.
```
sbatch -J  kely_hifiasm_default /blue/kawahara/yimingweng/universal_scripts/busco_default.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_hifisam_default/Kely_hifisam_default.fasta
```
- result: **C:96.6%[S:94.6%,D:2.0%],F:0.7%,M:2.7%,n:5286**

### **Purging Haplotigs**
1. Because the duplication is a bit too high, comparing to the published assemblies of other lepidopteran species which is usually around 1%, it is better to check and remove the duplicated haplotigs. Here I used [purge_haplotigs pipeline](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) to remove the duplicated haplotigs.
- step 1: map the raw reads to the genome assembly using [minimap2](https://github.com/lh3/minimap2). The script called `purge_minimap.slurm` is in the scripts folder. The output of this step includes the bam file from minimap2 and the histogram showing the read coverage for the contigs.
```
[yimingweng@login2 purging]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/purging

sbatch purge_minimap.slurm
```
<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/aligned.bam.histogram.png?raw=true?">

- step 2: I decide to apply cutoffs on read coverage, so determine the cutoffs of coverage for the contigs based on the coverage histogram. The cutoffs I would like to apply is min=5X and max=120X for contigs with 80% of reads falling out of this range. The script is called `purge_cutoff.slurm` in the scripts folder.
```
[yimingweng@login2 purging]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/purging

sbatch purge_cutoff.slurm
```

- step3: purge the haplotigs for the genome assembly with the information generated from step1 and step2.
```
[yimingweng@login2 purging]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/purging

sbatch purge_haplotigs.slurm
```
2. To check the assembly profile after purging steps, see how many contigs and nucleotides were removed.  
 -> contig count: **116 -> 62**  
 -> nucleotide count: **445.7Mbps -> 444Mbps**

3. Check the BUSCO result of updated genome assembly after purging haplotigs
```
sbatch -J  kely_hifiasm_purge /blue/kawahara/yimingweng/universal_scripts/busco_default.slurm /blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta
```
- result: **C:96.6%[S:95.5%,D:1.1%],F:0.7%,M:2.7%,n:5286**

### **Removing non-target sequences**
Although the assembly has gone through the haplotig purge pipeline, there is possibility that the assembly contains foreign (non-target) sequences. Blobplot from [blobtools](https://blobtools.readme.io/docs) is a good tool to detect and visualize the non-target sequence by examining the contig read depth, GC content, and blast results (to which taxon).
- step 1: blast the contigs to the nt (nucleotide) database from NCBI using megablast function in blastn tool. The script called `` is in the scripts folder. This step 
```
[yimingweng@login2 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

sbatch kely_megablast_nt.slurm
```
- step 2: prepare another input files for blobplot: the mapped reads file in bam format. Note that I have used minimap2 to do the similar work but that was on the original assembly. So I will have to repeat this work on the purged assembly again for the blobplot. The script called `kely_minimap_purge.slurm` was modified from `purge_minimap.slurm` to use the purged genome as input.
```
[yimingweng@login2 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

sbatch kely_minimap_purge.slurm
```
- step 3: run blobtools with the three input files:  
    (1) the genome assembly  
    (2) the reads mapped to the assembly in bam format   
    (3) the hit file by blasting the contig to the NCBI nt database
```
[yimingweng@login2 blobplot]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/blobplot

sbatch -J kely /blue/kawahara/yimingweng/universal_scripts/blobplot.slurm \
/blue/kawahara/yimingweng/Kely_genome_project/purging/Kely_purge_5X_120X.fasta \
/blue/kawahara/yimingweng/Kely_genome_project/blobplot/purging_aligned.bam \
/blue/kawahara/yimingweng/Kely_genome_project/blobplot/kely_purge_assembly.nt.mts1.hsp1.1e25.megablast.out \
kely_blobplot
```
<img src="https://github.com/yimingweng/Kely_genome_project/blob/main/blobplot/kely_blobplot.png?raw=true">

- there is a putative non-target contig in the assembly (ptg000079l), which has low read coverage and slightly higher GC content. Most importantly, it was blasted to Streptophyta, the plant, so remove it from the assembly. Now the updated assembly has:  
 -> contig count: **62 -> 61**  
 -> nucleotide count: **443,656,232bps -> 443,647,253bps**
- Since this contig doesn't contain any busco gene, so the busco result should remain same to the purged assembly.

### **Feature Annotation-1: repeats**
1. Usually the very first step of annotation is to identify the repeat regions using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler). Here I used the script called `repeatmodeler.slurm` to run repeatmodeler2. 
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

sbatch repeatmodeler.slurm /blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_final/kely_final_assembly.fasta kely_repeatmodeler
```

2. Use the output of repeatmodeler2 to run [RepeatMasker](https://github.com/rmhubley/RepeatMasker). Including the repeats from repeatmodeler2, I used 3 different evidence to mask the repeat regions for the genome:
    (1) mask the simple and short repeat, detected by RepeatMasker 
    (2) mask repeats based on existing databases (Repbase, Lepidoptera database)
    (3) mask genome based on the output of RepeatModeler2
Note that these scrip can be found in the scripts folder, and I used **soft-masking** for all the repeat regions.

- step1: simple and short repeat, script name: `kely_repeatmask_step1.slurm`
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

sbatch kely_repeatmask_step1.slurm
```

- step2: mask repeats based on existing databases (Repbase)[https://www.girinst.org/repbase/]
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

sbatch kely_repeatmask_step2.slurm
```
- step3: mask genome based on the output of RepeatModeler2
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation

sbatch kely_repeatmask_step3.slurm
```
- Calculate the masking percentage using the bash script called ``.
```
[yimingweng@login6 annotation]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/kely_repeatmasker_step34

bash /blue/kawahara/yimingweng/universal_scripts/maskrate.sh kely_final_assembly.fasta.masked.masked.masked
```

-> softmasking rate is **48.22%**  
-> hardmasking rate is 0%

### **Feature Annotation-2: gene model**
To get the best gene model, it is better to have RNAseq data so that the splicing sites and the gene features can be better caught by the predictor. However, I don't have RNAseq read for Keiferia, so the suboptimal alternative is to predict the gene model based on the ab initio (signal sensors and content sensors) and the protein sequence data from other insect species. 

1. In this case, I used [braker2](https://github.com/Gaius-Augustus/BRAKER#braker-with-proteins-of-any-evolutionary-distance) to get the gene model. There are 6 steps to finish this pipeline, and I have merged them into one script called `braker_prothint.slurm` in scripts folder.
- step1: download protein sequences from OrthoDB and take the repeat-masked genome as input, and run prohints to create protein database called prothint.gff
- step2: run Braker2 with the prohints database
- step3: use "selectSupportedSubsets.py" to extract the genes with full or partial support from the hints, and write it to a new gtf file called "supported_gene.gtf"
- step4: use supported_gene.gtf to generate a new protein fasta file called "braker_supported_gene.aa"
- step5: run busco without presence of isoform, further subset the braker_supported_gene.aa to a temporal fasta file (will be removed after the process finished) and this temporal fasta will be used to evaluate the completeness of the genome model with only the longest unique isoforms. That means the duplication from this busco result can not be referred to the presence of isoform.
- step6: get the final gft (kely_prohint_braker_final.gtf) from the anysupport.gff by fetching the gene sets with the list of supported transcripts 
```
[yimingweng@login6 kely_noRNA_model]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/prothint

sbatch -J kely_noRNA_model2 /blue/kawahara/yimingweng/universal_scripts/braker_prothint.slurm \
/blue/kawahara/yimingweng/Kely_genome_project/assemblies/kely_repeat_mask/kely_final_assembly_softmasked.fasta \
kely
```

- the result of this pipelines are the fasta file containing all the predicted gene sequences and the corresponding gff files. Generally, the files with name containing *anysupport* are the final products. Here is the busco result for the final gene model (after removing the smaller isoform):
    - Original BUSCO without removing unsupported gene (*not from this script*):   
    **C:96.2%[S:85.9%,D:10.3%],F:0.8%,M:3.0%,n:5286**
    - BUSCO result after removing the smaller isoforms:   
    **C:93.2%[S:91.7%,D:1.5%],F:0.7%,M:6.1%,n:5286**
    - The model comprises of **13841 genes** and **15405 transcripts**.


2. Run gFACs to get the annotation statistics
```
[yimingweng@login6 gfacs]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/gfacs

sbatch -J kely_final_model /blue/kawahara/yimingweng/universal_scripts/gfacs_statistics.slurm \
braker_2.1.2_gtf \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/prothint/kely_prohint_braker_final.gtf \
kely_gfacs \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/gfacs
```



### **Functional Annotation**
In this case, the gene model from the braker is called `kely_rerun2_anysupport_aa.fa`. I will blast those amino acid sequences in this file to different protein databases to get their predicted functions. For blasting job, I used [diamond](https://github.com/bbuchfink/diamond) to find the top 5 hits with the e-value cufoff to be 1e-5. The output file will be in the tsv format.

1. blast against the RefSeq non-redundant protein database (nr). I used a script called `diamond.slurm` to do run diamond. This script takes 4 arguments:
    (1) database for blast in fasta or dnmd format
    (2) protein sequences of the gene model from braker2 in fasta format
    (3) the e-value cutoff, usually 1*e-5
    (4) the output prefix for saving the tsv file
- download RefSeq non-redundant protein database (nr)
```
[yimingweng@login6 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/functional_annotation/diamond

wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
```
- run the diamond script
```
[yimingweng@login6 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/functional_annotation/diamond

sbatch -J kely_nr /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond/nr.dmnd /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa 0.00001 kely_nr_k5_1e5
```

2. blast against Swiss-Prot arthropod database (Reviewed Swiss-Prot for arthropod, data was downloaded through website, manually)
```
[yimingweng@login6 diamond]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/functional_annotation/diamond

sbatch -J kely_uniprot /blue/kawahara/yimingweng/universal_scripts/diamond.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/diamond/uniprot_arthropod.dmnd /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa 0.00001 kely_uniprot_k5_1e5
```

3. Usually we annotate functions for a gene model we want multiple lines of evidence. Here I am going to use [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/index.html) to annotate the genes. The script called `interproscan.slurm` is in the scripts folder.
```
[yimingweng@login5 interproscan]$ pwd
/blue/kawahara/yimingweng/Kely_genome_project/functional_annotation/interproscan

sbatch -J kely_interproscan /blue/kawahara/yimingweng/universal_scripts/interproscan.slurm /blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/kely_noRNA_model/kely_rerun2_anysupport_aa.fa kely_interproscan
```

4. additional annotation with [hmmer](http://hmmer.org/)
```
sbatch -J kely /blue/kawahara/yimingweng/universal_scripts/hmmer.slurm \
/blue/kawahara/yimingweng/Kely_genome_project/annotation/braker2/prothint/kely_anysupport_aa.fa \
kely
```

5. To annotate the protein sequences with GO terms, I used the [pannzer2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/).

6. The protein sequences can also be annotated with the gene function with the pathway it's involved. I use [KAAS](https://www.genome.jp/kegg/kaas/), the annotator for KEGG database. 

<br />

## References
Kawahara, A.Y., Storer, C.G. Markee, A., Heckenhauer, J., Powell, A., Plotkin, P., Hotaling, S., Cleland, T.P., Dikow, R.B., Dikow, T., Kuranishi, R.B., Messcher, R., Pauls, S.U., Stewart, R.J., Tojo, K., Frandsen, P.B. 2022. Long-read HiFi sequencing correctly assembles repetitive heavy fibroin silk genes in new moth and caddisfly genomes. Gigabyte, 2022  doi: 10.46471/gigabyte.64.

Li, X., Ellis, E., Plotkin, D., Imada, Y., Yago, M., Heckenhauer, J., Cleland, T., Dikow, R., Dikow, T., Storer, C.G., Kawahara, A.Y., Frandsen, P. 2021. First annotated genome of a mandibulate moth, Neomicropteryx cornuta, generated using PacBio HiFi sequencing. Genome Biology and Evolution 13(10). doi:10.1093/gbe/evab229. 

Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Marek Kokot, Maciej Długosz, Sebastian Deorowicz, KMC 3: counting and manipulating k-mer statistics, Bioinformatics, Volume 33, Issue 17, 01 September 2017, Pages 2759–2761, https://doi.org/10.1093/bioinformatics/btx304

Sebastian Deorowicz, Marek Kokot, Szymon Grabowski, Agnieszka Debudaj-Grabysz, KMC 2: fast and resource-frugal k-mer counting, Bioinformatics, Volume 31, Issue 10, 15 May 2015, Pages 1569–1576, https://doi.org/10.1093/bioinformatics/btv022

Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. Disk-based k-mer counting on a PC. BMC Bioinformatics 14, 160 (2013). https://doi.org/10.1186/1471-2105-14-160

Marek Kokot, Maciej Długosz, Sebastian Deorowicz, KMC 3: counting and manipulating k-mer statistics, Bioinformatics, Volume 33, Issue 17, 01 September 2017, Pages 2759–2761, https://doi.org/10.1093/bioinformatics/btx304

Sebastian Deorowicz, Marek Kokot, Szymon Grabowski, Agnieszka Debudaj-Grabysz, KMC 2: fast and resource-frugal k-mer counting, Bioinformatics, Volume 31, Issue 10, 15 May 2015, Pages 1569–1576, https://doi.org/10.1093/bioinformatics/btv022

Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. Disk-based k-mer counting on a PC. BMC Bioinformatics 14, 160 (2013). https://doi.org/10.1186/1471-2105-14-160

Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654

Manni, M., Berkeley, M. R., Seppey, M., & Zdobnov, E. M. (2021). BUSCO: Assessing genomic data quality and beyond. Current Protocols, 1, e323. doi: 10.1002/cpz1.323

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

Roach, M.J., Schmidt, S.A. & Borneman, A.R. Purge Haplotigs: allelic contig reassignment for third-gen diploid genome assemblies. BMC Bioinformatics 19, 460 (2018). https://doi.org/10.1186/s12859-018-2485-7

Laetsch DR and Blaxter ML. BlobTools: Interrogation of genome assemblies [version 1; peer review: 2 approved with reservations]. F1000Research 2017, 6:1287 (https://doi.org/10.12688/f1000research.12232.1)

Flynn, Jullien M., et al. "RepeatModeler2 for automated genomic discovery of transposable element families." Proceedings of the National Academy of Sciences 117.17 (2020): 9451-9457.

Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0.
2013-2015 <http://www.repeatmasker.org>.

Kim D. Pruitt, Tatiana Tatusova, Donna R. Maglott, NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins, Nucleic Acids Research, Volume 35, Issue suppl_1, 1 January 2007, Pages D61–D65, https://doi.org/10.1093/nar/gkl842

Wang Y, Wang Q, Huang H, Huang W, Chen Y, McGarvey PB, Wu CH, Arighi CN, UniProt Consortium.
A crowdsourcing open platform for literature curation in UniProt
Plos Biology. 19(12):e3001464 (2021)

Matthias Blum, Hsin-Yu Chang, Sara Chuguransky, Tiago Grego, Swaathi Kandasaamy, Alex Mitchell, Gift Nuka, Typhaine Paysan-Lafosse, Matloob Qureshi, Shriya Raj, Lorna Richardson, Gustavo A Salazar, Lowri Williams, Peer Bork, Alan Bridge, Julian Gough, Daniel H Haft, Ivica Letunic, Aron Marchler-Bauer, Huaiyu Mi, Darren A Natale, Marco Necci, Christine A Orengo, Arun P Pandurangan, Catherine Rivoire, Christian J A Sigrist, Ian Sillitoe, Narmada Thanki, Paul D Thomas, Silvio C E Tosatto, Cathy H Wu, Alex Bateman, Robert D Finn Nucleic Acids Research (2020), gkaa977, PMID: 33156333

Philip Jones, David Binns, Hsin-Yu Chang, Matthew Fraser, Weizhong Li, Craig McAnulla, Hamish McWilliam, John Maslen, Alex Mitchell, Gift Nuka, Sebastien Pesseat, Antony F. Quinn, Amaia Sangrador-Vegas, Maxim Scheremetjew, Siew-Yit Yong, Rodrigo Lopez, Sarah Hunter Bioinformatics (2014), PMID: 24451626

[HMMER](http://hmmer.org/): biosequence analysis using profile hidden Markov models

Törönen P, Medlar A, Holm L. PANNZER2: a rapid functional annotation web server. Nucleic Acids Res. 2018 Jul 2;46(W1):W84-W88. doi: 10.1093/nar/gky350. PMID: 29741643; PMCID: PMC6031051.

Moriya, Y., Itoh, M., Okuda, S., Yoshizawa, A., and Kanehisa, M.; KAAS: an automatic genome annotation and pathway reconstruction server. Nucleic Acids Res. 35, W182-W185 (2007).

Caballero M. and Wegrzyn J. 2019. gFACs: Gene Filtering, Analysis, and Conversion to Unify Genome Annotations Across Alignment and Gene Prediction Frameworks. Genomics_ Proteomics & Bioinformatics 17: 305-10