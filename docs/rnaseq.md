# The nf-core/rnaseq pipeline

In order to carry out a RNA-Seq analysis we will use the nf-core pipeline [rnaseq](https://nf-co.re/rnaseq/3.14.0).

## Overview

The pipeline is organised following the blocks we previously described: pre-processing, alignment (or pseudoalignment) and quantification and differential expression analysis.

![rnaseq overview](./img/nf-core-rnaseq_metro_map_grey.png)

In each process, the users can choose among a range of different options. Importantly, the users can decide to follow one of the two different routes in the alignment and quantification step:
- alignment and quantification (stage 2);
- pseudoalignment and quantification (stage 3). 

## Experimental Design


### Library design


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

In the following sections we will first prepare our references, then set our computational resources in order to be able to run the pipeline on a gitpod VM, edit the filtering settings and finally run the pipeline.

### Reference and annotation files

Following the considerations above, we will first of all edit the `nextflow.config` file in our working directory to add a new genome.
It is sufficient to add the following code to the `parameters` directive in the config.

```groovy
igenomes_base = '/workspace/gitpod/training/data/refs/'
genomes {
	'GRCh38chr21' {
        fasta                 = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta"
        fasta_fai             = "${params.igenomes_base}/sequence/Homo_sapiens_assembly38_chr21.fasta.fai"
        gff                   = "${params.igenomes_base}/trascriptome/gencode.v29.annotation_chr21_noversion.gff"
        transcript_fasta      = "${params.igenomes_base}/trascriptome/gencode.v29.transcripts_chr21_annotated.fa"
        salmon_index          = "${params.igenomes_base}/trascriptome/salmon_index_chr21.tar.gz"
	}
}
```

### Computing resources

Based on the choices we made when starting up the gitpod environment, we recommend to use the following additional parameters.
They can also be added to the parameters directive in the config file we just edited.

```groovy
params {
    max_cpus                  = 2
    max_memory                = '6.5GB'
    max_time                  = '2.h'
    use_annotation_cache_keys = true
}
```


### Filtering parameters



### Launching the pipeline

Now we are ready to launch the pipeline, and we can use the following command line:

```bash
nextflow run nf-core/rnaseq -r 3.12.0 \
--input  /workspace/gitpod/training//data/reads/rnaseq_samplesheet.csv \
--outdir . \
--genome GRCh38chr21 \
--pseudo_aligner salmon \
--skip_alignment \
--skip_biotype_qc \
-c rnaseq_nextflow.config \
-profile gls \
--skip_stringtie \
--skip_bigwig \
--skip_umi_extract \
--skip_trimming \
--skip_fastqc
```

