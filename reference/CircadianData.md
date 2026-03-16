# Circadian Liver Microarray Data (GSE11923)

This package provides processed liver transcriptome data from the Hughes
et al. (2009) study on circadian rhythms.

## Usage

``` r
data(eset)

data(eset_collapsed)

eset

eset_collapsed
```

## Format

Both objects are
[`ExpressionSet`](https://rdrr.io/pkg/Biobase/man/class.ExpressionSet.html)
objects:

- eset:

  Normalized, probe-level expression data (GCRMA). Contains 96 samples
  (48 time points in duplicate). Feature metadata includes gene symbols,
  gene names, and MAS5 presence/absence calls.

- eset_collapsed:

  Gene-level expression data. Probes were collapsed to unique gene
  symbols using the `MaxMean` method via WGCNA's `collapseRows`
  function.

An object of class `ExpressionSet` with 45101 rows and 48 columns.

An object of class `ExpressionSet` with 21032 rows and 48 columns.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11923>

## Details

The original data comes from GEO accession GSE11923. The data captures
the mouse liver transcriptome across a 48-hour period to study circadian
oscillations. The samples are labeled to reflect the ZT (Zeitgeber Time)
and replicate structure.

## References

Hughes, M. E., et al. (2009). Harmonics of Circadian Gene Transcription
in Mammals. *PLoS Genetics*, 5(4), e1000442.
