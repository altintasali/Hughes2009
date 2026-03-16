# Collapse probes to gene symbols

Collapses multiple probes mapping to the same gene symbol using the
WGCNA package.

## Usage

``` r
collapse_eset(eset, method = "MaxMean")
```

## Arguments

- eset:

  An ExpressionSet object.

- method:

  The WGCNA method to use. Options: 'MaxMean' (default), 'MaxVar',
  'Average', 'ME', 'First', 'Pattern'.

## Value

A new ExpressionSet object where rows are unique Gene Symbols.

## Examples

``` r
data(eset)
# Collapse probes to gene level using MaxMean method
gene_eset <- collapse_eset(eset, method = "MaxMean")
#> Collapsed 45101 probes into 21032 unique genes.
```
