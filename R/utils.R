#' @import Biobase
#' @import limma
#' @import WGCNA
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats setNames
NULL # This tells roxygen2 these imports belong to the package namespace

#' Filter ExpressionSet by MAS5 presence
#' 
#' Filters probes based on the number of samples where they were called 'Present' 
#' (MAS5 algorithm).
#'
#' @param eset An ExpressionSet object.
#' @param min_present Minimum number of samples (out of 48) where a probe 
#'   must be 'P' to be kept. Default is 10.
#' 
#' @return A filtered ExpressionSet object.
#' @export
#' @examples
#' data(eset)
#' # Filter to keep probes present in at least 12 samples
#' filtered_eset <- filter_mas5(eset, min_present = 12)
filter_mas5 <- function(eset, min_present = 10) {
  present_count <- rowSums(fData(eset)$MAS5_CALL == "P")
  eset[present_count >= min_present, ]
}

#' Collapse probes to gene symbols
#' 
#' Collapses multiple probes mapping to the same gene symbol using the WGCNA package.
#'
#' @param eset An ExpressionSet object.
#' @param method The WGCNA method to use. Options: 'MaxMean' (default), 
#'   'MaxVar', 'Average', 'ME', 'First', 'Pattern'.
#' 
#' @return A new ExpressionSet object where rows are unique Gene Symbols.
#' @export
#' @importFrom WGCNA collapseRows
#' @examples
#' data(eset)
#' # Collapse probes to gene level using MaxMean method
#' gene_eset <- collapse_eset(eset, method = "MaxMean")
collapse_eset <- function(eset, method = "MaxMean") {
  datExpr <- exprs(eset)
  gene_symbols <- fData(eset)$SYMBOL
  probe_ids <- featureNames(eset)
  mask <- !is.na(gene_symbols)
  
  # Perform collapse
  res <- WGCNA::collapseRows(datExpr[mask, ], 
                             gene_symbols[mask], 
                             probe_ids[mask], 
                             method = method)
  
  # Convert back to ExpressionSet
  new_eset <- ExpressionSet(assayData = res$datETcollapsed)
  pData(new_eset) <- pData(eset)
  
  message(sprintf("Collapsed %d probes into %d unique genes.", 
                  nrow(datExpr), nrow(res$datETcollapsed)))
  return(new_eset)
}

#' Plot Overlaid Density for all samples
#' @param eset An ExpressionSet object.
#' @param title Plot title.
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @import Biobase
plot_densities <- function(eset, title = "Sample Density Distributions") {
  dat <- exprs(eset)
  df <- reshape2::melt(dat)
  colnames(df) <- c("Probe", "Sample", "Intensity")
  
  ggplot(df, aes(x = .data$Intensity, group = .data$Sample)) +
    geom_density(alpha = 0.2) +
    theme_minimal() +
    labs(title = title, x = "Log2 Expression Intensity", y = "Density")
}

#' Plot Overlaid Boxplots for all samples
#' @param eset An ExpressionSet object.
#' @param title Plot title.
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @import Biobase
plot_boxplots <- function(eset, title = "Sample Expression Distribution") {
  dat <- exprs(eset)
  df <- reshape2::melt(dat)
  colnames(df) <- c("Probe", "Sample", "Intensity")
  
  ggplot(df, aes(x = .data$Sample, y = .data$Intensity)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = title, x = "Sample", y = "Log2 Expression Intensity")
}
