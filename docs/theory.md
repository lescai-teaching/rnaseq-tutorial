# RNAseq Analysis

Before we dive into the nf-core pipeline for analysing RNA-sequencing data, it's worth looking at some theoretical aspects of RNA-seq.

## Overview

Given the central role of RNA in a wide range of cellular and molecular functions, RNA-seq has emerged as a powerful tool for measuring the presence and levels of RNA species in biological samples. The technique is based on next-generation sequencing (NGS) technologies and is now considered the gold standard in the field of transcriptomics.

After RNA extraction and reverse transcription into complementary DNA (cDNA), the biological material is sequenced, generating NGS "reads" that correspond to the RNA captured in a specific cell, tissue, or organ at a given time. The sequencing data is then bioinformatically processed through a typical workflow summarised in the diagram below:

![overview](./img/Excalidraw_RNAseq.png)

In the scheme, we can identify three key phases in the workflow: data pre-processing, alignment and quantification, and differential expression analysis. In the data pre-processing step, the raw reads are processed to remove adapters and contaminants, and their quality is checked. Then, reads are mapped to a reference genome, and gene abundance is estimated. The workflow can also follow an alternative route based on lightweight alignment and quantification, reducing the time required for analysis. Finally, differentially expressed genes are identified using statistical tests, annotated, and visualised.

Depending on the user's needs, the workflow can include additional downstream analyses such as functional enrichment analysis, co-expression analysis, and integration with other omics data.

## Pre-processing 

The pre-processing of sequencing reads from RNA-seq data is a critical step to ensure the quality and accuracy of downstream analysis. The raw reads obtained from the sequencer are stored as [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files, which contain both the sequence data and quality scores. The initial processing step involves evaluating the quality of the reads and identifying potential issues such as adapter contamination, sequencing errors, or low-quality bases. The presence of adapters (short DNA sequences ligated to the ends of DNA fragments during library preparation) is detected through comparison with known adapter sequences or by using algorithms that identify adapter-like sequences, and removed in a process known as **read trimming**. Next, reads containing contaminants (genomic DNA and/or rRNA) and those with low-quality bases are filtered out. Finally, the quality of the filtered reads is checked again to ensure their suitability for downstream processing.

## Alignment (or lightweight alignment) and quantification

In the RNA-seq workflow, the alignment step refers to the process of mapping sequencing reads to a reference genome or transcriptome with the goal of determining the position and orientation of each read relative to the reference sequence.

Errors, gaps, or poor sequence quality regions, as well as insertions/deletions (INDELs), duplicated and repetitive regions in the reference sequence make this step challenging. Addressing these issues by choosing a high-quality reference and an appropriate aligner is essential for obtaining accurate results. A crucial component in the alignment step is the [annotation](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff) file, either in the form of General Feature Format file (GFF) or in form of Gene Transfer Format file (GTF). These files contain key information about the location and structure of genes and transcripts. For this reason they play a crucial role in both mapping sequencing reads accurately and in quantifying gene expression. Additionally, RNA-seq data often include reads that span exon-exon junctions, and the annotation files provide information about splice junctions allowing the inference of different isoforms.

The alignment and quantification steps can follow two different routes according to user preferences:
- alignment and quantification;
- lightweight alignment and quantification.

In the context of RNA-seq analysis, the alignment phase is often performed with splice-aware aligners that are utilised to take into account the splicing process. In addition to aligning reads across known splice junctions, splice-aware aligners also aim to detect novel splice junctions and alternative splicing events. Splice-aware aligners employ sophisticated algorithms and various optimization techniques, such as indexing the reference genome and parallel processing to achieve fast and scalable alignment. Popular splice-aware aligners include [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://github.com/DaehwanKimLab/hisat2). The following step is typically the quantification, which involves estimating the abundance (number of reads) assigned to each gene. Several tools are available to perform the quantification step, such as [featureCounts](https://subread.sourceforge.net/featureCounts.html), [HTSeq](https://htseq.readthedocs.io), [Salmon](http://combine-lab.github.io/salmon) and [RSEM](http://deweylab.github.io/RSEM).
The alignment and the quantification steps can be also performed with lightweight alignment tools, which include [Kallisto](https://pachterlab.github.io/kallisto/about.html), [Sailfish](http://www.cs.cmu.edu/~ckingsf/software/sailfish) and Salmon. These tools avoid a base-to-base alignment of reads providing quantification estimates faster than the classical splice-aware algorithms, but with a high accuracy. The resulting estimates are commonly referred to as **pseudocounts** or **abundance estimates**, which can be later utilised for downstream analysis.

## Differential expression (DE) analysis with DESeq2

The differential expression analysis is a statistical method to compare gene expression levels between different experimental conditions such disease vs. healthy (e.g., tumor tissue vs. healthy tissue), treatment vs control (e.g., sample treated with a specific stimulus, drug or compound vs untreated sample), tissue vs tissue (e.g., brain vs heart). The objective of differential expression analysis is to assess, for each gene, whether the observed differences in expression between groups are statistically significant, accounting for the variation observed within groups (replicates). This part of the analysis is typically performed in R with different packages that have been developed such as [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html). This tutorial focuses on **DESeq2**, a popular R package known for its robustness. For more detailed information and details refer to the [DESeq2 vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). The typical workflow is outlined in the flowchart below:

### Quality control

To ensure robust differential expression results, it’s common to start by exploring the sources of variation in the data. DESeq2 provides tools for quality control, including Principal Component Analysis (PCA) and hierarchical clustering. PCA is used to reduce the dimensionality of the data allowing visualisation of the samples in a lower-dimensional space. On the other hand, hierarchical clustering displays the correlation of gene expression for all pairwise combinations of samples in the dataset. These methods can identify sample groups with similar behavior, as well as potential outliers or unexpected patterns.
The quality control in DESeq2 typically uses **variance-stabilized (vst)** or **regularized log-transformed (rlog)** counts. Since raw counts in RNA-seq follow a discrete distribution, which is not suitable for many statistical and machine learning algorithms that assume continuous distributions, transformations like `vst` and `rlog` are applied to stabilise variance across genes. This stabilisation means that, after transformation, genes with both low and high expression levels will have variances that are more comparable, making the data more suitable for downstream analyses.
Finally, it’s often beneficial to filter out genes unlikely to show differential expression, such as those with zero counts across samples or extreme count outliers. By filtering out these genes, we can increase the sensitivity of our differential expression analysis and reduce the risk of false positive.

> [!NOTE]
> The `vst` or `rld` transformations are applied on normalised counts stored in the `dds` object, generated by running either the `DESeq()` or `estimateSizeFactors()` functions. Because size factors estimation is an early step of the `DESeq()` function, the transformations are applied immediately afterward.

### Design formula

The design formula specifies the sources of variation which the statistical needs to account for in modelling. It defines how the counts will be related to the experimental conditions, allowing the software to model the relationship between gene expression and the factors of interest, such as treatment, time points or batch effects. It is important to specify in the design formula the main factor of interest as well as additional covariates.
A well-defined design formula can account for potential confounding variables. For instance, if the experiment includes multiple conditions or batches, specifying these in the design helps to isolate the effect of the primary variable of interest on gene expression. An example is reported below:

```r
# Basic design with a single factor of interest

design = ~ condition

# Design including covariates to control for sex and developmental stage

design = ~ sex + developmental_stage + condition
```

In R, the tilde (`~`) is used in formula syntax to specify relationships between variables in statistical models. Here, it indicates that gene counts (dependent variable) will be modeled as a function of the variables specified (predictors).

> [!NOTE]
> The results will be not affected by the order of variables, but the common practice is to specificy the main source of variation in the last position of the design formula.


### Differential expression analysis

RNA-seq data typically include a large number of genes with low expression counts, reflecting that many genes are expressed at very low levels across samples. At the same time, RNA-seq data show a skewed distribution with a long right tail due to the absence of an upper limit for gene expression levels. This means that while most genes have low to moderate expression levels, a small number of genes are expressed at high levels. Accurate statistical modeling must therefore account for this distribution to avoid misleading conclusions.

   ![overview](./img/count_distribution.png)

The analysis starts with the input data generally composed by a matrix obtained in the alignment and quantification step that summarizes the expression levels of the different genes in each sample of the dataset. The rows of the matrix typically correspond to genes and the columns represent the samples. Another essential prerequisite is a metadata table describing the samples.
The core of the the differential expression analysis is the `DESeq()` function, a wrapper that streamlines multiple key steps in a single command. The different functions are reported here:

![overview](./img/DESeq_function.png)

> [!NOTE]
> While `DESeq()` combines these steps, a user could choose to carry out each function separately to have more control over the entire process.

The different steps are reported and explained in detail below:

1) **Normalisation**: since DESeq2 compares counts between sample groups for the same gene, it does not need to adjust for gene length. However, it is essential to account for variations in sequencing depth and RNA composition among samples. To normalise the data, DESeq2 utilizes size factors, which correct for these variations in library sizes and RNA composition.

The size factors are calculates using the **median ratio** method:

- **calculate the geometric mean**: for each gene, compute the geometric mean of its counts across all samples. This gives a row-wise geometric mean for each gene;

- **calculate ratios**: divide each gene's count by its geometric mean to obtain a ratio for each sample;

- **determine size factors**: for each sample, take the median of these ratios (column-wise) to derive the size factors;

- **normalise counts**: divide each raw count by the corresponding sample size factor to generate normalised counts.

The median ratio method assumes that not all genes are differentially expressed, so large outlier genes will not negatively influence the size factors. This approach is robust to imbalances in up-/down-regulation and can handle large numbers of differentially expressed genes.

> [!NOTE] 
> While normalised counts are useful for downstream visualisation of results, they should not be used as input for DESeq2. Instead, DESeq2 requires count data in the form of a matrix of integer values.

2) **Estimate dispersion and gene-wise dispersion**: the dispersion is a measure of how much the variance deviates from the mean. The estimation of the dispersion is essential to model the variance of the count data. Importatantly, RNA-seq data are characterised by **overdispersion**, where the variance in gene expression levels often exceeds the mean (variance > mean).

   ![overview](./img/overdispersion.png)

DESeq2 addresses this issue by employing the **negative binomial distribution**, which generalises the Poisson distribution by introducing an additional dispersion parameter. This parameter quantifies the extra variability present in RNA-seq data, providing a more realistic representation than the Poisson distribution, which assumes equal mean and variance.
DESeq2 starts by estimating the **common dispersion**, a single estimate of diepersion applicable to all genes on the dataset. This estimates provides a baseline for variability across all genes in the dataset. Next, DESeq2 estimates **gene-wise dispersion**, separate estimate of dispersion for each individual gene, taking into account that different genes may exhibit varying levels of expression variability due to biological differences. The dispersion parameter (α) is related to the mean (μ) and variance of the data, as described by the equation:

`Var = μ + α ⋅ μ2`

A key feature of DESeq2's dispersion estimates is their negative correlation with the mean and positive correlation with variance. So, genes with low expression have higher dispersion values and genes with high expression lower tend to have lower dispersion. Additionally, genes sharing similar mean expression levels can exhibit different dispersion estimates based on their variance. To improve the accuracy of dispersion estimates, DESeq2 assumes that genes with similar expression profiles share similar dispersion patterns and leverages this information to refine the estimates.

- **Mean-dispersion relationship**: this process, known as dispersion fitting, models the relationship between the mean expression level of a gene and its dispersion. In this process, DESeq2 identifies a trend in the dispersion estimates across genes. The fitted curve, typically a smooth curve, describes how dispersion changes as a function of the mean expression level.

- **Final dispersion estimates**: DESeq2 refines the gene-wise dispersion by shrinking it toward the fitted curve. The "shrinkage" helps control for overfitting, particularly in genes with low counts or few replicates, and makes the dispersion estimates more reliable.However, genes with exceptionally high dispersion values are not shrinked because they probably deviate from the modelling assumption, exhibiting elevated variability due to biological or technical factors. Shrinking these values could lead to false positives.

![overview](./img/DESeq_dispersion_estimates.png)

**Fitting model and testing**: DESeq2 fits a generalized linear model (GLM) to the normalised counts using the calculated size factors and final dispersion estimates.
The GLM fit is then used to perform hypothesis testing for differential expression using either a Wald test or a likelihood ratio test. However, when performing multiple tests, such as in the case of RNA-seq data where thousands of genes are tested, the risk of false positives increases. To account for this, DESeq2 employs multiple test correction methods (Benjamini-Hochberg procedure is the default) to adjust the p-values and control the false discovery rate (FDR). 


> [!NOTE] 
> For example by setting the FDR cutoff to < 0.05, 5% of genes identified as differentially expressed are expected to be false positive. For instance, if 400 genes are identified as differentially expressed with an FDR cutoff of 0.05, you would expect 20 of them to be false positives.


After identifying DE genes using DESeq2, it is essential to interpret the biological significance of these genes through functional analysis. This involves examining the roles of the differentially expressed genes in various biological processes, molecular functions and pathways providing insights into the underlying mechanisms driving the observed changes in gene expression. Different tools are available to carry out these functional analyses such as [Gene Ontology](https://geneontology.org), [Reactome](https://reactome.org/), [KEGG](https://www.genome.jp/kegg), [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), [g:Profiler](https://biit.cs.ut.ee/gprofiler), and [WikiPathways](https://www.wikipathways.org).
