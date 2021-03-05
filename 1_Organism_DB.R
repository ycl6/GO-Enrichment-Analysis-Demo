#' ---
#' title: "GO Enrichment Analysis (Part 1: About Organism Annotation DB)"
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
#' **Bioconductor AnnotationData Packages:** http://www.bioconductor.org/packages/release/data/annotation/
#' 
#' **AnnotationHub:**: https://bioconductor.org/packages/AnnotationHub/
#' 
#' **License:** GPL-3.0
#' 
#' # Introduction
#' 
#' There are many organism-level (`org`) packages readily available on Bioconductor. They provide mappings between a central identifier (e.g. Entrez Gene identifiers) and other identifiers (e.g. ensembl ID, Refseq Identifiers, GO Identifiers, etc). 
#' 
#' The name of an `org` package is always of the form `org.<Sp>.<id>.db` (e.g. `org.Hs.eg.db`) where `<Sp>` is a 2-letter abbreviation of the organism (e.g. Hs for *Homo sapiens*) and `<id>` is an abbreviation (in lower-case) describing the type of central identifier (e.g. `eg` for gene identifiers assigned by the Entrez Gene, or `sgd` for Saccharomyces Genome Database). Most of the Bioconductor annotation packages are updated every 6 months.
#' 
#' # Start `R`
#' 
#' ```
#' cd /ngs/GO-Enrichment-Analysis-Demo
#' 
#' R
#' ```
#' 
#' # Using `BiocManager`
#' 
#' List available organism-level packages for installation in BiocManager.
#' 
## ----available-org------------------------------------------------------------
BiocManager::available("^org\\.")

#' 
#' ## Install *Arabidopsis* `org` package
#' 
#' As an example, let's download and install the *Arabidopsis thaliana* (thale cress) package.
#' 
## ----install-db-BiocManager---------------------------------------------------
BiocManager::install("org.At.tair.db")
library(org.At.tair.db)
org.At.tair.db

#' 
#' # Using `AnnotationHub`
#' 
#' Above method returns a limited number of organism-level annotation packages. There are a lot more packages available from the Bioconductor's [AnnotationHub](https://bioconductor.org/packages/AnnotationHub/) service.
#' 
#' To search, download and install packages from the AnnotationHub service, install `AnnotationHub` if it is not yet installed in your machine.
#' 
## ----install-AnnotationHub----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("AnnotationHub")

#' 
#' ## Create an `AnnotationHub` object
#' 
## ----load-AnnotationHub-------------------------------------------------------
library(AnnotationHub)

ah <- AnnotationHub()

# URL for the online AnnotationHub
hubUrl(ah)

#' 
#' ## Summary of available records
#' 
## ----print-summary------------------------------------------------------------
ah

# Number of resources
length(ah)

#' 
#' ## Query the hub for `org` records
#' 
#' Search for organism-level packages with a pattern-matching string "`^org\\.`".
#' 
## ----query-AnnotationHub------------------------------------------------------
db <- query(ah, "^org\\.")
df <- mcols(db)
class(df)

#' 
#' ## Show query results
#' 
#' Show query results stored in `DFrame`.
#' 
## ----query-df-----------------------------------------------------------------
# Column names
cbind(colnames(df))

# Number of org records
nrow(df)

# Show df
df[,c("title", "species")]

#' 
#' ## Download *Felis* `org` package
#' 
#' Let's search and install the *Felis catus* (cat) package.
#' 
## ----install-db-AnnotationHub-------------------------------------------------
# Search df with keyword
data.frame(df[grep("Felis", df$species), c("title", "species", "rdatadateadded")])

# Retrieve package with for "Felis catus"
rn <- rownames(df[df$species == "Felis catus",])
org.Fc.eg.db <- ah[[rn]]
org.Fc.eg.db

#' 
#' ## Show record status
#' 
## ----show-record-status-------------------------------------------------------
recordStatus(ah, rn)

#' 
#' ## Load from local cache
#' 
#' After retrieving an annotation package, it will be placed in the local AnnotationHub cache. You can used it again without having to download the package. 
#' 
## -----------------------------------------------------------------------------
# Location of the local AnnotationHub cache
hubCache(ah)

# Load from cache
org.Fc.eg.db <- ah[[rn]]

#' 
#' ## Clear local cache
#' 
#' You can use the `removeCache` function to removes all local AnnotationHub database and all related resources. 
#' 
## ----remove-cache, eval = FALSE-----------------------------------------------
## removeCache(ah, ask = TRUE)

#' 
#' # Discover `org` db objects
#' 
#' ## `columns`
#' 
#' Shows which kinds of data can be returned for the `AnnotationDb` object. 
#' 
#' Both objects contain **Gene Ontology** mapping information.
#' 
## ----discover-columns---------------------------------------------------------
columns(org.At.tair.db)
columns(org.Fc.eg.db)

#' 
#' ## `keytypes`
#' 
#' Shows which columns can be used as keys.
#' 
## ----discover-keytypes--------------------------------------------------------
keytypes(org.At.tair.db)
keytypes(org.Fc.eg.db)

#' 
#' ## `keys`
#' 
#' Returns values (or keys) that can be expected for a given keytype. By default it will return the primary keys for the database.
#' 
## ----discover-keys------------------------------------------------------------
head(keys(org.At.tair.db), 10)	# Primary keys
head(keys(org.At.tair.db, keytype = "SYMBOL"), 10)
head(keys(org.At.tair.db, keytype = "GO"), 10)

head(keys(org.Fc.eg.db), 10)	# Primary keys
head(keys(org.Fc.eg.db, keytype = "SYMBOL"), 10)
head(keys(org.Fc.eg.db, keytype = "GO"), 10)

#' 
#' ## `select`
#' 
#' Retrieve the data as a `data.frame` based on parameters for selected `keys`, `columns` and `keytype` arguments.
#' 
#' ### Ex1: Given TAIR ID, retrieves SYMBOL
#' 
## ----discover-select-1--------------------------------------------------------
myKeys <- head(keys(org.At.tair.db, keytype = "TAIR"), 10)
myKeys
select(org.At.tair.db, keys = myKeys, columns = "SYMBOL", keytype = "TAIR")

#' 
#' ### Ex2: Given SYMBOL, retrieves ENTREZID ID
#' 
## ----discover-select-2--------------------------------------------------------
myKeys <- c("CCA1", "LHY", "PRR7", "PRR9") # morning loop components
select(org.At.tair.db, keys = myKeys, columns = "ENTREZID", keytype = "SYMBOL")

#' 
#' ### Ex3: Given ENSEMBL ID, retrieves SYMBOL
#' 
## ----discover-select-3--------------------------------------------------------
myKeys <- head(keys(org.Fc.eg.db, keytype = "ENSEMBL"), 10)
myKeys
select(org.Fc.eg.db, keys = myKeys, columns = "SYMBOL", keytype = "ENSEMBL")

#' 
#' ### Ex4: Given SYMBOL, retrieves ENSEMBL ID and ENTREZID ID
#' 
## ----discover-select-4--------------------------------------------------------
myKeys <- c("ASIP", "MC1R") # coat color patterns
select(org.Fc.eg.db, keys = myKeys, columns = c("ENSEMBL", "ENTREZID"), keytype = "SYMBOL")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

#' 
