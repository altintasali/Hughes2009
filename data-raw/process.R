#==============================================================================#
# Project info ---- 
#==============================================================================#
# Project    : Hughes 2009
# Data       : Microarray data
# Task       : Pre-process 
# Author     : Ali Altintas 
# Date       : 16.03.2026
#==============================================================================#

rm(list = ls(all = TRUE)) #Remove Workspace

#==============================================================================#
# Libraries ----
#==============================================================================#
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load("usethis", "GEOquery", "gcrma", "affy", "mouse4302.db", 
#                "Biobase", "WGCNA")
library(GEOquery)
library(gcrma)
library(affy)
library(gcrma)
library(mouse4302.db)
library(Biobase)
library(WGCNA)

#==============================================================================#
# Pre-processing ----
#==============================================================================#
# 1. Download and Extract
getGEOSuppFiles("GSE11923")
untar("GSE11923/GSE11923_RAW.tar", exdir = "GSE11923/CEL")
cels <- list.files("GSE11923/CEL", pattern = "CEL", full.names = TRUE)

# 2. Read and immediate cleanup of raw files
raw_data <- ReadAffy(filenames = cels)
unlink("GSE11923", recursive = TRUE)

# 3. Normalization & MAS5 Calls
eset <- gcrma(raw_data)
mas5_calls <- mas5calls(raw_data)
present_matrix <- exprs(mas5_calls)

# 4. Metadata Download and immediate cleanup
gse_meta_list <- getGEO("GSE11923", destdir = ".")
eset_meta <- gse_meta_list[[1]]
pData(eset) <- pData(eset_meta)

# Cleanup metadata file(s)
file.remove(list.files(pattern = "^GSE11923.*\\.txt\\.gz$"))
file.remove(list.files(pattern = "^GPL1261.*\\.soft\\.gz$"))

# 5. Annotate Probes to Gene Symbols
probe_ids <- featureNames(eset)
anno <- select(mouse4302.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
anno_unique <- anno[!duplicated(anno$PROBEID), ]
rownames(anno_unique) <- anno_unique$PROBEID
fData(eset) <- anno_unique[probe_ids, ]
fData(eset)$MAS5_CALL <- present_matrix[probe_ids, ]

# 6. Create Collapsed version (MaxMean method)
datExpr <- exprs(eset)
gene_symbols <- fData(eset)$SYMBOL
mapping_mask <- !is.na(gene_symbols)

collapsed_all <- collapseRows(datExpr[mapping_mask, ], 
                              gene_symbols[mapping_mask], 
                              probe_ids[mapping_mask], 
                              method = "MaxMean")

eset_collapsed <- ExpressionSet(assayData = collapsed_all$datETcollapsed)
pData(eset_collapsed) <- pData(eset)

# 7. Save only the 2 required objects
usethis::use_data(eset, eset_collapsed, overwrite = TRUE, compress = "xz")
