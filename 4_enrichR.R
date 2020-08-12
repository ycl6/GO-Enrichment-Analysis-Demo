#' ---
#' title: "GO Enrichment Analysis (Part 4: Enrichr & enrichR)"
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
#' **Enrichr - a gene set enrichment analysis web server:** [Website](https://amp.pharm.mssm.edu/Enrichr/), [Paper](https://doi.org/10.1093/nar/gkw377)
#' 
#' **enrichR - an R interface to all 'Enrichr' databases:** [CRAN](https://CRAN.R-project.org/package=enrichR)
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
library("enrichR")
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
up.genes <- data[up.idx,]$GeneSymbol
dn.genes <- data[dn.idx,]$GeneSymbol

#' 
## -----------------------------------------------------------------------------
head(up.genes, 10)
head(dn.genes, 10)

#' 
#' Alternatively, if you only have Ensembl gene ID
#' 
## -----------------------------------------------------------------------------
up.genes <- data[up.idx,]$GeneID
dn.genes <- data[dn.idx,]$GeneID

#' 
## -----------------------------------------------------------------------------
head(up.genes, 10)
head(dn.genes, 10)

#' 
#' # Prepare input data
#' 
#' We would need to convert any other identifier format to `SYMBOL` which is the required input identifier format. This can be done by using the `select` function from `AnnotationDbi` that we saw in [Part 1](1_Organism_DB.html) of this demo, or by using the "Biological Id TRanslator" `bitr` function from `clusterProfiler` which is a wrapper function of `AnnotationDbi::select`. 
#' 
#' Here, we will use `bitr` here to see how this can be done.
#' 
## ----use-bitr-----------------------------------------------------------------
# Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
up.genes.df = clusterProfiler::bitr(up.genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
head(up.genes.df, 10)

dn.genes.df = clusterProfiler::bitr(dn.genes, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Mm.eg.db")
head(dn.genes.df, 10)

#' 
#' # Connecting to Enrichr web service
#' 
#' List available databases from Enrichr
#' 
## ----list-dbs-----------------------------------------------------------------
dbs <- listEnrichrDbs()
dbs <- dbs[order(dbs$libraryName),]

class(dbs)
dim(dbs)
head(dbs)

#' 
#' Show all database names.
#' 
## ----show-dbs-----------------------------------------------------------------
dbs$libraryName

#' 
#' Search for mouse databases with keyword `"Mouse"`
#' 
## ----show-mouse-dbs-----------------------------------------------------------
dbs[grep("Mouse",dbs$libraryName),]$libraryName

#' 
#' # Select databases
#' 
## ----select-dbs---------------------------------------------------------------
dbs_go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
dbs_pw <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "BioPlanet_2019")
dbs_dd <- c("PheWeb_2019", "ClinVar_2019")

#' 
#' # Perform enrichment analysis
#' 
#' ## GO analysis
#' 
## ----run-enrichr-go-----------------------------------------------------------
upEnriched_go <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_go)
dnEnriched_go <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_go)

#' 
## -----------------------------------------------------------------------------
class(upEnriched_go)
names(upEnriched_go)

head(upEnriched_go[[1]], 10)

#' 
#' ## Pathway analysis
#' 
## ----run-enrichr-pw-----------------------------------------------------------
upEnriched_pw <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_pw)
dnEnriched_pw <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_pw)

#' 
## -----------------------------------------------------------------------------
class(upEnriched_pw)
names(upEnriched_pw)

head(upEnriched_pw[[1]], 10)

#' 
#' ## Diseases/Drugs analysis
#' 
## ----run-enrichr-dd-----------------------------------------------------------
upEnriched_dd <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_dd)
dnEnriched_dd <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_dd)

#' 
## -----------------------------------------------------------------------------
class(upEnriched_dd)
names(upEnriched_dd)

head(upEnriched_dd[[1]], 10)

#' 
#' # Create functions
#' 
#' ## `.enrichment_prep_df` function
#' 
## ----enrichment_prep_df-------------------------------------------------------
.enrichment_prep_df <- function(df, showTerms, orderBy) {

    if(is.null(showTerms)) {
	showTerms = nrow(df)
    } else if(!is.numeric(showTerms)) {
        stop(paste0("showTerms '", showTerms, "' is invalid."))
    }

    Annotated <- as.numeric(sub("^\\d+/", "", as.character(df$Overlap)))
    Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))

    # Build data frame
    df <- cbind(df, data.frame(Annotated = Annotated, Significant = Significant, 
			       stringsAsFactors = FALSE))

    # Order data frame (P.value or Combined.Score)
    if(orderBy == "Combined.Score") {
        idx <- order(df$Combined.Score, decreasing = TRUE)
    } else {
        idx <- order(df$P.value, decreasing = FALSE)
    }
    df <- df[idx,]

    # Subset to selected number of terms
    if(showTerms <= nrow(df)) {
        df <- df[1:showTerms,]
    }

    return(df)
}

#' 
#' ## `plotEnrich` function
#' 
#' Here, we will use a `plotEnrich` function to visualise Enrichr results as bar plots.
#' 
## ----plotEnrich---------------------------------------------------------------
plotEnrich <- function(df, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value", 
		       xlab = NULL, ylab = NULL, title = NULL) {

    if(!is.numeric(numChar)) {
        stop(paste0("numChar '", numChar, "' is invalid."))
    }

    df <- .enrichment_prep_df(df, showTerms, orderBy)

    # Create trimmed name (as seen in topGO)
    shortName <- paste(substr(df$Term, 1, numChar),
                       ifelse(nchar(df$Term) > numChar, '...', ''), sep = '')
    df$shortName = shortName
    df$shortName <- factor(df$shortName, levels = rev(unique(df$shortName)))
    df$Ratio <- df$Significant/df$Annotated

    # Define fill variable (P.value or Combined.Score)
    if(orderBy == "Combined.Score") {
	fill <- "Combined.Score"
    } else {
       	fill <- "P.value"
    }

    # Define y variable (Count or Ratio)
    if(y != "Ratio") {
        y <- "Significant"
    }

    # Define variable mapping
    map <- aes_string(x = "shortName", y = y, fill = fill)

    # Define labels
    if(is.null(xlab)) {
        xlab <- "Enriched terms"
    }

    if(is.null(ylab)) {
        if(y == "Ratio") {
            ylab <- "Gene ratio"
        } else {
            ylab <- "Gene count"
        }
    }

    if(is.null(title)) {
        title <- "Enrichment analysis by Enrichr"
    }

    # Make the ggplot
    p <- ggplot(df, map) + geom_bar(stat = "identity") + coord_flip() + theme_bw()

    if(orderBy == "Combined.Score") {
        p <- p + scale_fill_continuous(low = "blue", high = "red") + 
		guides(fill = guide_colorbar(title = "Combined Score", reverse = FALSE))
    } else {
	p <- p + scale_fill_continuous(low = "red", high = "blue") + 
		guides(fill = guide_colorbar(title = "P value", reverse = TRUE))
    }

    # Adjust theme components
    p <- p + theme(axis.text.x = element_text(colour = "black", vjust = 1), 
		   axis.text.y = element_text(colour = "black", hjust = 1),
		   axis.title = element_text(color = "black", margin = margin(10, 5, 0, 0)),
		   axis.title.y = element_text(angle = 90))

    p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)

    return(p)
}

#' 
#' ## `printEnrich` function
#' 
#' Here, we will use a `printEnrich` function to output Enrichr results to text files.
#' 
## ----printEnrich--------------------------------------------------------------
printEnrich <- function(data, prefix = "enrichr", showTerms = NULL, columns = c(1:9)) {

    if(!is.numeric(columns)) {
        stop(paste0("columns '", columns, "' is invalid."))
    }

    for (i in 1:length(data)) {
        dbname <- names(data)[i]
        df <- data[[i]]

	df <- .enrichment_prep_df(df, showTerms, orderBy = "P.value")
	df <- df[, !colnames(df) %in% c("Annotated", "Significant")]

        if(any(columns > ncol(df))) {
            stop("Undefined columns selected")
        }

	filename <- paste0(prefix, "_", dbname, ".txt")
	write.table(df, file = filename, sep = "\t", quote = F, row.names = F, col.names = T)
    }
}

#' 
#' # Plot enrichment
#' 
#' Demonstrate using different paramters to plot enrichment.
#' 
## ----plot-results, fig.width = 8, fig.height = 6, fig.align = "center", dpi = 100----
plotEnrich(upEnriched_go[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
plotEnrich(upEnriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value")
plotEnrich(upEnriched_dd[[2]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score")

#' 
#' # Output results to files
#' 
## ----output-results-----------------------------------------------------------
printEnrich(upEnriched_go, prefix = "enrichr-GO-up", showTerms = 20)
printEnrich(dnEnriched_go, prefix = "enrichr-GO-dn", showTerms = 20)
printEnrich(upEnriched_pw, prefix = "enrichr-PW-up", showTerms = 20)
printEnrich(dnEnriched_pw, prefix = "enrichr-PW-dn", showTerms = 20)
printEnrich(upEnriched_dd, prefix = "enrichr-DD-up", showTerms = 20)
printEnrich(dnEnriched_dd, prefix = "enrichr-DD-dn", showTerms = 20)

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

