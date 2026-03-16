# Filter ExpressionSet by MAS5 presence

Filters probes based on the number of samples where they were called
'Present' (MAS5 algorithm).

## Usage

``` r
filter_mas5(eset, min_present = 10)
```

## Arguments

- eset:

  An ExpressionSet object.

- min_present:

  Minimum number of samples (out of 48) where a probe must be 'P' to be
  kept. Default is 10.

## Value

A filtered ExpressionSet object.

## Examples

``` r
data(eset)
# Filter to keep probes present in at least 12 samples
filtered_eset <- filter_mas5(eset, min_present = 12)
```
