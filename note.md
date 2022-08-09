### 08/05/2020
- Getting to know that we will have PacBio sequence reads for Tomato Pinworm, *Keiferia lycopersicella* and Groosefoot Groundling Moth, *Scrobipalpa atriplicella* back in a week or two. 

- The computational work will be done through the cluster machine in UF called HiperGator where an account is required.

- It might be a good idea to get familiar with this species and the genome assembly pipeline for the PacBio reads.   

<br />

### 08/08/2020
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

