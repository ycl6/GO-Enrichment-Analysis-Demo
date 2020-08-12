# GO Enrichment Analysis Demo

### Demonstration of Gene Ontology (GO) enrichment analysis

-----

**GO-Enrichment-Analysis-Demo web version:** https://ycl6.github.io/GO-Enrichment-Analysis-Demo

**topGO:** [Bioconductor](https://bioconductor.org/packages/topGO/), [Paper](https://doi.org/10.1093/bioinformatics/btl140)

**topGO (with `enrichment_barplot` function):** [GitHub](https://github.com/ycl6/topGO-feat)

**clusterProfiler:** [Bioconductor](https://bioconductor.org/packages/clusterProfiler/), [Paper](https://doi.org/10.1089/omi.2011.0118), [Documentation](https://yulab-smu.github.io/clusterProfiler-book/)

**enrichplot:** [Bioconductor](https://bioconductor.org/packages/enrichplot/)

**Enrichr:** [Website](https://amp.pharm.mssm.edu/Enrichr/), [Paper](https://doi.org/10.1093/nar/gkw377)

**enrichR:** [CRAN](https://CRAN.R-project.org/package=enrichR)

**Demo Dataset:** [E-MTAB-8411](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411) from *[The clock gene Bmal1 inhibits macrophage motility, phagocytosis, and impairs defense against pneumonia.         ](https://doi.org/10.1073/pnas.1915932117) PNAS. 2020;117(3):1543-1551.*

**License:** GPL-3.0

## Introduction

In this tutorial, I will use the sequencing data from [E-MTAB-8411](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411) to demonstrate how to use `topGO` and `clusterProfiler` to perform GO enrichment   analysis. 

This RNA-seq dataset contains 5 wild-type (GK1, GK3, GK5, GK7 and GK9) and 4 macrophage-specific *Bmal1* knockout samples (GK2, GK4, GK6 and GK10). The fastq files were retreived from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411/samples/) and processed with the follow steps:

1. Adapter and QC trimming with `BBDuk`
2. Mapping to Mouse genome (assemble MM10 and [GENCODE M24](https://www.gencodegenes.org/mouse/release_M24.html) annotation) with `STAR`
3. Differential expression analysis with `DESeq2`

We will perform the enrichment analysis using output from `DESeq2`.

## Prerequisites

### Demo Data

Download the result table from GitHub

```S
cd /ngs/GO-Enrichment-Analysis-Demo

wget https://raw.githubusercontent.com/ycl6/GO-Enrichment-Analysis-Demo/master/DESeq2_DEG.txt
```

### R Packages

```R
# From CRAN
install.packages(c("devtools", "data.table"))

# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("org.Mm.eg.db", "clusterProfiler", "enrichplot"))

# From GitHub
# topGO with enrichment_barplot()
devtools::install_github("ycl6/topGO-feat", ref = "v2.41.0-barplot")
```

## Resources

### 1. About Organism Annotation DB (1_Organism_DB.R; [HTML](https://ycl6.github.io/GO-Enrichment-Analysis-Demo/1_Organism_DB.html))

### 2. topGO (2_topGO.R; [HTML](https://ycl6.github.io/GO-Enrichment-Analysis-Demo/2_topGO.html))

### 3. clusterProfiler (3_clusterProfiler.R; [HTML](https://ycl6.github.io/GO-Enrichment-Analysis-Demo/3_clusterProfiler.html))

### 4. Enrichr & enrichR (4_enrichR.R; [HTML](https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html))

