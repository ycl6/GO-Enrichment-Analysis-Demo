#' ---
#' title: "GO Enrichment Analysis (Part 3: clusterProfiler)"
#' description: |
#'   Demonstration of Gene Ontology (GO) enrichment analysis
#' author:
#'   - name: "I-Hsuan Lin"
#'     url: https://github.com/ycl6
#'     affiliation: University of Manchester
#'     affiliation_url: https://www.manchester.ac.uk/
#' date: '`r format(Sys.Date(), "%B %d, %Y")`'
#' output:
#'     rmarkdown::html_document:
#'         theme: united
#'         highlight: tango
#'         self_contained: true
#'         toc: true
#'         toc_depth: 2
#'         toc_float:
#'             collapsed: false
#'             smooth_scroll: true
#' ---
#' 

#' 
#' -----
#' 
#' **clusterProfiler:** [Bioconductor](https://bioconductor.org/packages/clusterProfiler/), [Paper](https://doi.org/10.1089/omi.2011.0118), [Documentation](https://yulab-smu.github.io/clusterProfiler-book/)
#' 
#' **enrichplot:** [Bioconductor](https://bioconductor.org/packages/enrichplot/)
#' 
#' **Demo Dataset:** [E-MTAB-8411](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411) from *[The clock gene Bmal1 inhibits macrophage motility, phagocytosis, and impairs defense against pneumonia.         ](https://doi.org/10.1073/pnas.1915932117) PNAS. 2020;117(3):1543-1551.*
#' 
#' **License:** GPL-3.0
#' 
#' 
#' # Start `R`
#' 
#' ```
#' cd /ngs/GO-Enrichment-Analysis-Demo
#' 
#' R
#' ```
#' 
#' # Load package and set path
#' 
## ----load-libraries, message = FALSE------------------------------------------
library("clusterProfiler")
library("enrichplot")
library("org.Mm.eg.db")
library("ggplot2")

#' 
#' # Load data
#' 
#' If you have downloaded the `DESeq2_DEG.txt` file with `wget`:
#' 
## ----load-data----------------------------------------------------------------
data <- data.table::fread("DESeq2_DEG.txt")
data$GeneID <- substr(data$GeneID, 1, 18)

#' 
#' If you like to donwload the file in `R` now:
#' 
## ----eval = FALSE-------------------------------------------------------------
## data <- data.table::fread("https://raw.githubusercontent.com/ycl6/GO-Enrichment-Analysis-Demo/master/DESeq2_DEG.txt")
## data$GeneID <- substr(data$GeneID, 1, 18)

#' 
## -----------------------------------------------------------------------------
data

#' 
#' # Define significance threshold
#' 
## ----significance-threshold---------------------------------------------------
up.idx <- which(data$padj < 0.05 & data$log2fc > 0)	# FDR < 0.05 and logFC > 0
dn.idx <- which(data$padj < 0.05 & data$log2fc < 0)	# FDR < 0.05 and logFC < 0

#' 
## -----------------------------------------------------------------------------
dim(data)
length(up.idx)
length(dn.idx)

#' 
#' # Define significant genes
#' 
## ----significant-genes--------------------------------------------------------
all.genes <- data$GeneSymbol
up.genes <- data[up.idx,]$GeneSymbol
dn.genes <- data[dn.idx,]$GeneSymbol

#' 
## -----------------------------------------------------------------------------
head(up.genes, 10)
head(dn.genes, 10)

#' 
#' Alternatively, if you only have Ensembl gene ID
#' 
## ----eval = FALSE-------------------------------------------------------------
## all.genes <- data$GeneID
## up.genes <- data[up.idx,]$GeneID
## dn.genes <- data[dn.idx,]$GeneID

#' 
#' # Decide the sub-ontology to test
#' 
#' - BP: Biological Process
#' - CC: Cellular Component
#' - MF: Molecular Function
#' 
## ----define-ontology----------------------------------------------------------
ontology <- "BP"

#' 
#' # Set outfile prefix
#' 
## ----outfile-prefix-----------------------------------------------------------
outTitle <- paste0("clusterProfiler_GO-", ontology, "_ORA_simplify")
outTitle

#' 
#' # Prepare input data
#' 
#' We use all the annotated genes (`all.genes`) as the "gene universe" and the differentially expressed genes (`up.genes` and `dn.genes`) as our genes of interest. 
#' 
#' We would need to convert any other identifier format to `ENTREZID` which is the required input identifier format. This can be done by using the `select` function from `AnnotationDbi` that we saw in [Part 1](1_Organism_DB.html) of this demo, or by using the "Biological Id TRanslator" `bitr` function from `clusterProfiler` which is a wrapper function of `AnnotationDbi::select`. 
#' 
#' Here, we will use `bitr` here to see how this can be done.
#' 
## ----use-bitr-----------------------------------------------------------------
# Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
all.genes.df = bitr(all.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(all.genes.df, 10)

up.genes.df = bitr(up.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(up.genes.df, 10)

dn.genes.df = bitr(dn.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
head(dn.genes.df, 10)

#' 
#' # Perform enrichment analysis
#' 
#' - P-value cutoff on enrichment tests is set at 0.05.
#' - P-values adjustment method is set to "BH" (Benjamini & Hochberg, same as "fdr").
#' - Q-value cutoff is set at 0.05. Results that must pass (1) `pvalueCutoff` on unadjusted p-values, (2) `pvalueCutoff` on adjusted p-values and (3) `qvalueCutoff` on q-values.
#' 
## ----run-enrichGO-------------------------------------------------------------
upEGO = enrichGO(gene = up.genes.df$ENTREZID, universe = all.genes.df$ENTREZID, 
		 OrgDb = "org.Mm.eg.db", ont = ontology, pvalueCutoff = 0.05, 
		 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
upEGO

dnEGO = enrichGO(gene = dn.genes.df$ENTREZID, universe = all.genes.df$ENTREZID,
                 OrgDb = "org.Mm.eg.db", ont = ontology, pvalueCutoff = 0.05,
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
dnEGO

#' 
#' # Reduce GO term redundancy 
#' 
#' We then use the `simplify` function to reduce redundancy of enriched GO terms. We will used the default parameters to run the function.
#' 
## ----run-simplify-------------------------------------------------------------
upSimGO = simplify(upEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
		   semData = NULL)
nrow(upSimGO)

dnSimGO = simplify(dnEGO, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
		   semData = NULL)
nrow(dnSimGO)

#' 
#' # Plot enrichment
#' 
#' We will use the `barplot` function from `enrichplot` (previously a function in `clusterProfiler`)
#' 
#' ### Up-regulated genes
#' 
## ----plot-enrichment-up-------------------------------------------------------
png(paste0(outTitle, "_up.png"), width = 9, height = 6, units = "in", res = 300)
barplot(upSimGO, showCategory = 20) + 
	ggtitle(paste0("GO-", ontology," ORA of up-regulated genes")) + 
	xlab("Enriched terms") + ylab("Count")
invisible(dev.off())

#' 

#' 
#' ### Down-regulated genes
#' 
## ----plot-enrichment-down-----------------------------------------------------
png(paste0(outTitle, "_dn.png"), width = 9, height = 6, units = "in", res = 300)
barplot(dnSimGO, showCategory = 20) + 
	ggtitle(paste0("GO-", ontology," ORA of down-regulated genes")) +
        xlab("Enriched terms") + ylab("Count")
invisible(dev.off())

#' 

#' 
#' # Create result table
#' 
#' Save `enrichGO` results as a `data.frame`
#' 
## ----result-table-------------------------------------------------------------
up.tab = upSimGO@result
dn.tab = dnSimGO@result

#' 
#' # Output results to files
#' 
## ----output-results-----------------------------------------------------------
write.table(up.tab, file = paste0(outTitle, "_up.txt"), sep = "\t", quote = F, 
	    row.names = F, col.names = T)

write.table(dn.tab, file = paste0(outTitle, "_dn.txt"), sep = "\t", quote = F, 
	    row.names = F, col.names = T)

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

