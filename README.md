
# Hughes2009: Circadian Liver Transcriptome

[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://altintasali.github.io/Hughes2009/)

The `Hughes2009` package provides processed liver transcriptome data
from the Hughes et al. (2009) study on circadian rhythms. It contains
normalized expression data and gene-level summaries, designed for
testing circadian analysis workflows and co-expression network
construction.

# Installation

You can install the development version of `Hughes2009` from GitHub:

``` r
install.packages("devtools")
devtools::install_github("altintasali/Hughes2009")
```

# Quick Start

This package provides an `ExpressionSet` object containing mouse liver
microarray data across a 48-hour period.

``` r
library(Hughes2009)
library(Biobase)
```

``` r
# Load the dataset
data(eset)

# Peek at the data
show(eset)
```

    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 45101 features, 48 samples 
    ##   element names: exprs 
    ## protocolData
    ##   sampleNames: GSM301348.CEL.gz GSM301349.CEL.gz ... GSM301395.CEL.gz
    ##     (48 total)
    ##   varLabels: ScanDate
    ##   varMetadata: labelDescription
    ## phenoData
    ##   sampleNames: GSM301348 GSM301349 ... GSM301395 (48 total)
    ##   varLabels: title geo_accession ... relation.1 (33 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1415670_at 1415671_at ... AFFX-r2-P1-cre-5_at (45101
    ##     total)
    ##   fvarLabels: PROBEID SYMBOL GENENAME MAS5_CALL
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation: mouse4302

# Vignettes

``` r
# View the vignette
browseVignettes("Hughes2009")
```

# Full Documentation

[Read the full documentation and tutorials
online](https://altintasali.github.io/Hughes2009/)
