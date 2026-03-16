#' Circadian Liver Microarray Data (GSE11923)
#'
#' This package provides processed liver transcriptome data from the 
#' Hughes et al. (2009) study on circadian rhythms. 
#'
#' @docType data
#' @format Both objects are \code{\link[Biobase]{ExpressionSet}} objects:
#' \describe{
#'   \item{eset}{Normalized, probe-level expression data (GCRMA). Contains 96 
#'   samples (48 time points in duplicate). Feature metadata includes gene 
#'   symbols, gene names, and MAS5 presence/absence calls.}
#'   \item{eset_collapsed}{Gene-level expression data. Probes were collapsed to 
#'   unique gene symbols using the \code{MaxMean} method via WGCNA's 
#'   \code{collapseRows} function.}
#' }
#' 
#' @details 
#' The original data comes from GEO accession GSE11923. The data captures 
#' the mouse liver transcriptome across a 48-hour period to study circadian 
#' oscillations. The samples are labeled to reflect the ZT (Zeitgeber Time) 
#' and replicate structure.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11923}
#' @references Hughes, M. E., et al. (2009). Harmonics of Circadian Gene 
#' Transcription in Mammals. *PLoS Genetics*, 5(4), e1000442.
#' 
#' @usage data(eset)
#' @usage data(eset_collapsed)
#' @name CircadianData
NULL

#' @rdname CircadianData
"eset"

#' @rdname CircadianData
"eset_collapsed"