# The nf-core/rnaseq pipeline

In order to carry out a RNA-Seq analysis we will use the nf-core pipeline [rnaseq](https://nf-co.re/rnaseq/3.14.0).

## Overview

The pipeline is organised following the blocks we previously described: pre-processing, alignment (or pseudoalignment) and quantification and differential expression analysis.

![rnaseq overview](./img/nf-core-rnaseq_metro_map_grey.png)

In each process, the users can choose among a range of different options. Importantly, the users can decide to follow one of the two different routes in the alignment and quantification step:
- alignment and quantification (stage 2);
- pseudoalignment and quantification (stage 3). 

## Experimental Design

In RNA-seq experimental design, the number of reads and the number of biological replicates are two critical factors that need to be carefully considered. While it may seem intuitive that having a large number of reads is always desirable, this is not necessarily the case. In fact, having an excessive number of reads can lead to unnecessary costs and computational burdens, without providing significant improvements in data quality or biological insight. Instead, it is often more beneficial to prioritize the number of biological replicates, as these allow for the capture of natural biological variation and provide a more accurate representation of the underlying biological phenomenon. Biological replicates involve collecting and sequencing RNA from separate biological samples from different individuals, different tissues, or different time points. They are essential for identifying genuine changes in gene expression and distinguishing them from technical noise. In contrast, technical replicates involve sequencing the same RNA sample multiple times to assess the technical variability of the sequencing platform and ensure consistency in the results.
To ensure the validity of the results, it is crucial to carefully consider the trade-off between the number of biological replicates and sequencing depth. Higher sequencing depth can improve the ability to detect lowly expressed genes and provide more precise quantification of gene expression. However, beyond a certain threshold, increasing sequencing depth yields diminishing returns in terms of data quality and biological insight, making it important to find an optimal balance. Therefore, it is often more beneficial to prioritize biological replicates over simply increasing sequencing depth or technical replicates.
Furthermore, the number of reads and biological replicates should be estimated through a calculation of the statistical power, which takes into account the desired fold change, alpha level, and effect size. This calculation can inform the experimental design and provide a basis for setting the log2 fold change (log2FC) threshold in the subsequent differential expression (DE) analysis. By incorporating multiple biological replicates into the experimental design, optimizing sequencing depth, and estimating the required number of reads and replicates through statistical power calculations, it is possible tp increase the statistical power of the analysis, reduce the risk of false positives, and gain a more comprehensive understanding of the biological system being studied.

### Library design

RNA-seq library design involves a range of critical decisions, including the choice between paired-end and single-end sequencing strategies. This choice is influenced by several key factors, such as fragment size selection, strand specificity preservation, and optimal read length determination. Paired-end sequencing, where both ends of a fragment are sequenced, provides valuable information about structural variations and transcript isoforms and can improve mapping accuracy, especially for longer transcripts and repetitive regions. In contrast, single-end sequencing, where only one end of the fragment is sequenced, can be more cost-effective while still providing high-quality data for gene expression analysis. The decision between paired-end and single-end sequencing ultimately depends on the research question and experimental goals. For instance, paired-end sequencing may be preferred if the focus is on identifying novel transcripts or characterizing transcript isoforms. However, if the primary goal is to quantify gene expression levels, single-end sequencing may be sufficient. The type of RNA being sequenced also plays a crucial role in library design, as different RNA species (e.g., mRNA, total RNA, small RNA) require distinct protocols. Additionally, the chosen read length can impact the ability to detect specific features, such as splice junctions or repeated regions. Therefore, the choice between paired-end and single-end reads in RNA-seq library design depends on the specific research goals, the nature of the RNA samples, budget constraints, alongside considerations for data analysis complexity and computational resources.

### Reference genome

nf-core pipelines make use of the Illumina iGenomes collection as [reference genomes](https://nf-co.re/docs/usage/reference_genomes).
Before starting the analysis, the users might want to check whether the genome they need is part of this collection.
They also might want to consider downloading the reference locally, when running on premises: this would be useful for multiple runs and to speed up the analysis. In this case the parameter `--igenomes_base` might be used to pass the root directory of the downloaded references. 

One might also need to use custom files: in this case the users might either provide specific parameters at command line, or create a config file adding a new section to the `genome` object. See [here](https://nf-co.re/docs/usage/reference_genomes#custom-genomes) for more details.

We will follow this specific approach in this tutorial, since the data we will be using have been simulated on chromosome 21 of the Human GRCh38 reference, and we have prepared fasta, indexes and annotation files containing only this chromosome locally.

### Input files

The input data should be provided in a CSV file, according to a format that is largely common for nf-core pipelines.
The format is described in the [rnaseq usage page](https://nf-co.re/rnaseq/3.14.0/docs/usage).

## Running nf-core/rnaseq

In the following sections we will:
- prepare our references;
- set our computational resources in order to be able to run the pipeline on a gitpod VM;
- edit the filtering settings;
- run the pipeline.

To prevent losing connection, open the **Terminal** and type:

```bash
sleep 2h
```

This command will keep the gitpod session active for exactly 2 hours, providing sufficient time to complete our analysis without interruption.

![overview](./img/sleep_2h.png)

Now come back to our **Console**.

### Reference and annotation files

Following the considerations above, we will first of all edit the `nextflow.config` file in our working directory to add a new genome.
It is sufficient to add the following code to the `parameters` directive in the config.

```groovy
igenomes_base = '/workspace/gitpod/training/data/refs/'
genomes {
        'GRCh38chr21' {
            fasta                 = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta"
            fasta_fai             = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta.fai"
            gff                   = "${params.igenomes_base}/gencode_v29_chr21_parsed.gff"
            transcript_fasta      = "${params.igenomes_base}/gencode.v29.transcripts_chr21.fa"
            salmon_index          = "${params.igenomes_base}/salmon_index_chr21.tar.gz"
	}
}
```

### Computing resources

Based on the choices we made when starting up the gitpod environment, we recommend to use the following additional parameters.
They can also be added to the parameters directive in the config file we just edited.

```groovy
params {
    max_cpus      = 2
    max_memory    = '6.5GB'
    max_time      = '2.h'
}
```

### Launching the pipeline

Now we are ready to launch the pipeline, and we can use the following command line:

```bash
nextflow run nf-core/rnaseq -r 3.12.0 \
--input /workspace/gitpod/training/data/reads/rnaseq_samplesheet.csv \
--outdir ./results_star_salmon \
--genome GRCh38chr21 \
--aligner star_salmon \
--pseudo_aligner salmon \
--skip_biotype_qc \
--skip_stringtie \
--skip_bigwig \
--skip_umi_extract \
--skip_trimming \
--skip_fastqc \
--skip_markduplicates \
--skip_dupradar \
--skip_rseqc \
--skip_qualimap
```
