*TBXT* chromatin looping analysis in lung cancer
================
<h4>
Authors: <i>Reuben M. Yaa</i>
</h4>
<h4>
README updated: <i>Aug-22-2022</i>
</h4>

## Introduction

While high levels of TBXT have been detected in metastatic lung cancer,
We used circularised chromosome capture combined with sequencing
(4C-seq) to determine and physical genome wide *cis* -interactions
differences around *TBXT* loci in 1) lung cancer cells (H460), a TBXT
expressing cell line 2) lung cancer cells (H358 and A549), non-TBXT
expressors

## 4C template and library

4C libraries were prepared in replicates using *Nla*III and *Dpn*II
restriction enzymes and sequenced and sequenced at 150bp PE

## 4C data analysis

Mapping

Reads were processed using [4C
pipeline](https://github.com/mangiruben/TBXT_regulation/blob/main/code/4C_pipeline.pl)

Interaction landscape callings

Statistically interacting fragments were called using [interact\_calling
pipeline](https://github.com/mangiruben/TBXT_regulation/blob/main/code/interact_calling.R)

## Downstream analysis

Analysis of the statistically interacting fragments

The fragments now called loops were analysed in different ways as
follows;

-   Venn diagrams for comparing loops between *TBXT* expression
    conditions
-   Merged loop for visualisation of the regulation landscape linear
    genome using
    [spider\_plots](https://github.com/mangiruben/TBXT_regulation/blob/main/code/spider_plot.R)
-   Clustering of loops as shared, gained or lost using Bedtools suit
-   Density analysis of clustered loops
-   Motif analysis of gained loops

Epigenome enrichment analysis

ChIP-seq analysis of H358, A549 and H460 publicly available datasets
using
[ChIP\_pipeline](https://github.com/mangiruben/TBXT_regulation/blob/main/code/chip_seq_SE.pl)

<hr>
