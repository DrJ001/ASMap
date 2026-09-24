# Individual marker and interval statistics for an R/qtl cross object

Individual marker and interval statistics for an R/qtl cross object

## Usage

``` r
statMark(cross, chr, stat.type = c("marker","interval"),
         map.function = "kosambi")
```

## Arguments

- cross:

  An qtl `cross` object with class structure `"bc"`, `"dh"`, `"riself"`,
  `"bcsft"`. (see
  [`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  for more details.)

- chr:

  Character vector of linkage group names used for subsetting the
  linkage map.

- stat.type:

  Character string of either `"marker"` or `"interval"` or both.
  `"marker"` produces individual marker related statistics and
  `"interval"` produces interval related statistics for the current map
  order (see Details).

- map.function:

  Character string of either `"kosambi"`, `"haldane"`, `"morgan"` or
  `"cf"` defining the map function to be used for interval related
  statistics.

## Details

If `"marker"` is chosen then a call to `geno.table` from qtl is used to
return individual marker statistics for segregation distortion, as well
as allele and missing value proportions. For the current map order the
number of double crossovers at each marker are also returned.

If `"interval"` is chosen then interval statistics are returned for the
current map order. These include the estimated recombination fraction
and LOD score between adjacent markers, calculated from `est.rf` in qtl.
Also returned are the map interval distances and converted map
recombination fractions extracted from the `"map"` component of each
linkage group as well as the actual number of recombinations between
markers.

This function is used in `profileMark` to plot any combination of
returned linkage map statistics on a single graphical display.

## Value

A list named by the `stat.type` used in the call. Each element is a data
frame of statistics with columns named by the statistic.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`profileMark`](https://drj001.github.io/ASMap/reference/profileMark.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## produce all statistics

sm <- statMark(mapDH, stat.type = c("marker","interval"))
```
