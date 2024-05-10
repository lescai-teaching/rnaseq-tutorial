# The nf-core/rnaseq pipeline

In order to carry out a RNA-Seq analysis we will use the nf-core pipeline [rnaseq](https://nf-co.re/rnaseq/3.14.0).

## Overview

The pipeline is organised following the blocks we previously described: pre-processing, alignment (or pseudoalignment) and quantification and differential expression analysis.

![rnaseq overview](./img/nf-core-rnaseq_metro_map_grey.png)

In each process, the user can choose among a range of different options. Importantly, the user can decide to follow one of the two different routes in the alignment and quantification step:
- alignment and quantification (stage 2);
- pseudoalignment and quantification (stage 3). 

## Experimental Design

nf-core pipelines make use of the Illumina iGenomes collection as [reference genomes](https://nf-co.re/docs/usage/reference_genomes).


### Library design


### Reference genome


### Input files

The input data should be provided in a CSV file, according to a format that is largely common for nf-core pipelines.
The format is described in the [rnaseq usage page](https://nf-co.re/rnaseq/3.14.0/docs/usage).

## Running nf-core/rnaseq


### Reference and annotation files



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

```

