# Interpretation

Once the differential expression genes have been identified, the next crucial step is to interpret the results. This involves examining the tables and plots generated during the analysis to understand the biological implications of the data. In this part of the tutorial, we will delve into the results by discussing the significant genes identified, their expression patterns, and the visual representations of these findings. We will explore various tables and plots, such as volcano plots, MA plots, and heatmaps, to gain insights into the underlying biological processes and validate the reliability of our results.

The first plot we will examine is the Principal Component Analysis (PCA) plot. Since we're working with simulated data, our metadata is relatively simple, consisting of just three variables: sample, condition, and replica. In a typical RNA-seq experiment, however, metadata can be complex and encompass a wide range of variables that could contribute to sample variation, such as sex, age and developmental stage. 

![overview](./img/PCA.png)

By plotting the PCA on the PC1 and PC2 axes, using [condition] as the main variable of interest, we can quickly identify the primary source of variation in our data. By accounting for this variation in our design model, we should be able to detect more differentially expressed genes related to [condition]. When working with real data, it's often useful to plot the data using different variables to explore how much variation is explained by the first two PCs. Depending on the results, it may be informative to examine variation on additional PC axes, such as PC3 and PC4, to gain a more comprehensive understanding of the data.

Next, we will examine the hierarchical clustering plot to explore the relationships between samples based on their gene expression profiles. The heatmap is organized such that samples with similar expression profiles are close to each other, allowing us to identify patterns and structures in the data.

![overview](./img/hierarchical_clustering.png)

We can observe a high degree of correlation in the plot. Remember that to create this plot, we utilized the [dist()] function, so in the legend on the right, a value of 0 corresponds to high correlation, while a value of 5 corresponds to very low correlation. Similar to PCA, we can see that samples tend to cluster together according to [condition].

Overall, the integration of these plots suggests that we are working with high-quality data and can confidently proceed to the differential expression analysis.

From this point, we will examine plots that are generated after running the differential expression analysis. These plots are not quality control (QC) plots, but rather plots that help us interpret the results of the differential expression analysis. 
After running the [results()] function, a good way to start to have an idea about the results of our differential expression analysis is to look at the MA plot. 

![overview](./img/MA_plot.png)

By default, genes are coloured in blue if the padj is less than 0.1 and the log2foldchange greater than or less than 0. Genes that fall outside the plotting region are represented as open triangles. Note that we have not yet applied a filter to select only significant, which we define as those with a padj value less than 0.5 and a log2 fold change of at least 1 or -1.

After filtering our genes of interest according to our threshold, let's have a look to our significatnt genes

```bash
log2 fold change (MLE): condition treatment vs control 
Wald test p-value: condition treatment vs control 
DataFrame with 5 rows and 6 columns
gene            baseMean      log2FoldChange    lfcSE       stat       pvalue          padj
ENSG00000205726 121653.2815   2.895764          0.1501028   19.291872   6.285623e-83    1.477121e-80
ENSG00000142192  51025.6955   3.023104          0.1868371   16.180420   6.931957e-59    8.145049e-57
ENSG00000142156  20789.8570   2.975533          0.2127267   13.987585   1.856093e-44    1.453940e-42
ENSG00000159231    458.8502   -1.195471         0.3170778   -3.770276   1.630670e-04    9.580189e-03
ENSG00000156282    481.6481   1.091949          0.3084051   3.540632    3.991703e-04    1.876100e-02
```

After the identification of DE genes, it's informative to visualize the expression of specific genes of interest. Using the [plotCounts()] function directly on the [dds] object allows us to examine individual gene expression profiles without accessing the full [res] object.

![overview](./img/plotCounts.png)

In our example, post-treatment, we observe a significant increase in the expression of the ENSG00000142192 gene, highlighting its responsiveness to the experimental conditions.

To gain a comprehensive overview of the transcriptional profile, the volcano plot represents a highly informative tool.

![overview](./img/volcanoplot.png)

The treatment induced differential expression in five genes, with one downregulated and four upregulated. This plot visually represents the numerical results reported in the table above.

Finally, we can create a heatmap using the normalized expression counts of genes identified as significantly differentially expressed. The resulting heatmap visualizes how the expression of significant genes varies across samples. Each row represents a gene, and each column represents a sample. The color intensity in the heatmap reflects the normalized expression levels: red colors indicate higher expression, while blue colors indicate lower expression.

![overview](./img/heatmap_de_genes.png)

By examining the heatmap, we can visually identify the expression patterns of our five significant differentially expressed genes. This visualization allows us to identify how these genes respond to the treatment. The heatmap provides a clear and intuitive way to explore gene expression dynamics and can guide further investigation into the *DESEq* processes underlying observed changes.

## Over Representation Analysis (ORA)

Finally, we can attempt to assign biological significance to our differentially expressed genes through Enrichment Analysis (ORA). The ORA analysis identifies specific biological pathways, molecular functions, and cellular processes that are enriched with our differentially expressed genes.

![overview](./img/enrichment_plot.png)

Based on the results of the ORA, we identified two differentially expressed genes that appear to be involved in redox reduction, specifically in the reduction of carbonyl groups. Furthermore, these genes seem to play a role in the endoplasmic reticulum during clathrin-mediated endocytosis and endosome recycling. This suggests that they may be involved in metabolic pathways related to detoxification and/or xenobiotic metabolism.

Note that this analysis is based on simulated data and serves only as a demonstration to illustrate the classical workflow of the analysis.

## Conclusions

In this tutorial, we have walked through the steps of performing differential expression analysis using DESeq2, from preparing the data to interpreting the results. We have seen how to identify differentially expressed genes, visualize the results, and perform enrichment analysis to gain insights into the biological significance of the results. By following this tutorial, you should now be able to perform differential expression analysis using DESeq2 and interpret the results in the context of your own research question. Remember to always carefully evaluate the quality of your data and the assumptions of the analysis, and to consider the biological relevance of the results in the context of your research.