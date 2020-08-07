#' ---
#' title: "GO Enrichment Analysis (Part 2: topGO)"
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
#' **topGO:** [Bioconductor](https://bioconductor.org/packages/topGO/), [Paper](https://doi.org/10.1093/bioinformatics/btl140)
#' 
#' **topGO (with `enrichment_barplot` function):** [GitHub](https://github.com/ycl6/topGO-feat)
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
library("topGO")
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
#' # Decide test algorithm
#' 
#' The default algorithm used by the topGO package is `weight01`, it is a mixture between the `elim` and the `weight` algorithms. Possible choices includes:
#' 
#' - classic
#' - elim
#' - weight
#' - weight01
#' - lea
#' - parentchild
#' 
## ----define-algorithm---------------------------------------------------------
algorithm <- "weight01"

#' 
#' # Define the statistical test used
#' 
#' For tests based on gene counts: 
#' 
#' - Fischer's exact test (`statistic = "fisher"`)
#' 
#' For tests based on gene scores or gene ranks: 
#' 
#' - Kolmogorov-Smirnov test (`statistic = "ks"`)
#' - t-test (`statistic = "t"`)
#' 
#' For tests based on gene expression:
#' 
#' - globaltest (`statistic = "globaltest"`)
#' 
#' 
## ----define-statistic---------------------------------------------------------
statistic <- "fisher"

#' 
#' # Set outfile prefix
#' 
## ----outfile-prefix-----------------------------------------------------------
outTitle <- paste0("topGO_GO-", ontology, "_ORA_", algorithm,"_", statistic)
outTitle

#' 
#' # Prepare input data
#' 
#' We use all the annotated genes (`all.genes`) as the "gene universe" or `geneList` and the differentially expressed genes (`up.genes` and `dn.genes`) as our genes of interest. So, we will define the `geneList` as 1 if a gene is differentially expressed, 0 otherwise
#' 
#' 
## ----prepare-input-up---------------------------------------------------------
upList <- factor(as.integer(all.genes %in% up.genes))
names(upList) <- all.genes

head(upList, 30)
table(upList)

#' 
## ----prepare-input-dn---------------------------------------------------------
dnList <- factor(as.integer(all.genes %in% dn.genes))
names(dnList) <- all.genes

head(dnList, 30)
table(dnList)

#' 
#' # Create topGOdata object
#' 
#' Here, our `geneList` is named with gene symbol, hence we use `ID = "SYMBOL"`. If `geneList` is named with Ensembl gene ID, you need to use `ID = "ENSEMBL"`
#' 
## ----create-topGOdata---------------------------------------------------------
upGOdata <- new("topGOdata", ontology = ontology, allGenes = upList,geneSel = function(x)(x == 1), 
		nodeSize = 10, annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "SYMBOL")

dnGOdata <- new("topGOdata", ontology = ontology, allGenes = dnList,geneSel = function(x)(x == 1), 
		nodeSize = 10, annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "SYMBOL")

#' 
#' # Test for enrichment
#' 
## ----test-enrichment----------------------------------------------------------
upRes <- runTest(upGOdata, algorithm = algorithm, statistic = statistic)
upRes

dnRes <- runTest(dnGOdata, algorithm = algorithm, statistic = statistic)
dnRes

#' 
#' # Plot enrichment
#' 
#' ### Up-regulated genes
#' 
## ----plot-enrichment-up-------------------------------------------------------
png(paste0(outTitle, "_up.png"), width = 8, height = 6, units = "in", res = 300)
enrichment_barplot(upGOdata, upRes, showTerms = 20, numChar = 50, orderBy = "Scores", 
		   title = paste0("GO-", ontology," ORA of up-regulated genes"))
invisible(dev.off())

#' 

#' 
#' ### Down-regulated genes
#' 
## ----plot-enrichment-down-----------------------------------------------------
png(paste0(outTitle, "_dn.png"), width = 8, height = 6, units = "in", res = 300)
enrichment_barplot(dnGOdata, dnRes, showTerms = 20, numChar = 50, orderBy = "Scores",
                   title = paste0("GO-", ontology," ORA of down-regulated genes"))
invisible(dev.off())

#' 

#' 
#' # Create result table
#' 
#' Here, we retrieve test results of top 20 GO terms, you can use `topNodes = length(usedGO(GOdata_Obj))` to retrieve results of all available GO terms
#' 
## ----GenTable-----------------------------------------------------------------
up.tab <- GenTable(upGOdata, Pval = upRes, topNodes = 20)
up.tab

dn.tab <- GenTable(dnGOdata, Pval = dnRes, topNodes = 20)
dn.tab

# Update table with full GO term name
up.tab$Term <- sapply(up.tab$"GO.ID", function(go) Term(GO.db::GOTERM[[go]]))
up.tab$Term

dn.tab$Term <- sapply(dn.tab$"GO.ID", function(go) Term(GO.db::GOTERM[[go]]))
dn.tab$Term

#' 
#' # Add gene symbol to result table
#' 
#' The `GenTable` function only provide the significant gene count for each significant GO term, here we will add the gene symbols to the result table
#' 
## ----add-gene-symbols---------------------------------------------------------
# Obtain the list of significant genes
up.sigGenes <- sigGenes(upGOdata)
dn.sigGenes <- sigGenes(dnGOdata)

# Retrieve gene symbols for each GO from the test result
up.AnnoList <- lapply(up.tab$"GO.ID", 
		      function(x) as.character(unlist(genesInTerm(object = upGOdata, whichGO = x))))
dn.AnnoList <- lapply(dn.tab$"GO.ID", 
		      function(x) as.character(unlist(genesInTerm(object = dnGOdata, whichGO = x))))

up.SigList <- lapply(up.AnnoList, function(x) intersect(x, up.sigGenes))
dn.SigList <- lapply(dn.AnnoList, function(x) intersect(x, dn.sigGenes))

# Coerce gene list to a comma-separated vector
up.tab$Genes <- sapply(up.SigList, paste, collapse = ",")
dn.tab$Genes <- sapply(dn.SigList, paste, collapse = ",")

#' 
## -----------------------------------------------------------------------------
cbind(head(up.tab$Genes, 5))
cbind(head(dn.tab$Genes, 5))

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

