---
title: "TBXT chromatin loopoing analysis in lung cancer"  
author: "<h4>Authors: <i>Reuben M.Yaa, Rafal D. Acamel, Juan Tena</i></h4>" 
date: "<h4>README updated: <i>`r format( Sys.Date(), '%b-%d-%Y')`</i></h4>"
output:
  github_document
---

## Introduction
While high levels of TBXT have been detected in metastatic lung cancer, 
We used circularised chromosome capture combined with sequencing (4C-seq) to determine physical  genome wide cis -interactions from TBXT promoter in 
1)lung cancer cells H460, a *TBXT* expressing cell line compared to 
2)H358 and A549, non-*TBXT* expressors  


## 4C template and library
4C libraries were prepared in replicates using NlaIII and DpnII restriction enzymes and sequenced read length 150bp PE

*4C data mapping*  
Reads were processed using 4C pipeline  

*Interaction lanscape callings*  
Statisticlly interacting fragments were called using interact_calling.R
 
## Downstream analysis*  
Difference analysis of the statistically interacting fragments(now called loops) was done
 - Venn diagrams for comparing loops bettween conditions
 - Merged loop for visualisation of the regulation landscape linear genome using spider plot
 - Clustering of loops as shared, gained or lost using Bedtools suit
 - Density analysis of clustered loops
 - Motif analysis of gained loops
 
 *Epigenome enrichment analysiss*  
 ChIP-seq analysis of H358, A549 and H460 publicly available datasets
<hr>
