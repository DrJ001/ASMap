# Profile individual genotype statistics for an R/qtl cross object

Profile individual genotype statistics for the current linkage map order
of and R/qtl cross object

## Usage

``` r
profileGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo",
           "miss"), id = "Genotype", xo.lambda = NULL, ...)
```

## Arguments

- cross:

  An qtl `cross` object with class structure
  `"bc", "dh", "riself", "bcsft"`. (see
  [`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  for more details.)

- chr:

  Character vector of linkage group names used for subsetting the
  linkage map.

- bychr:

  Logical vector determining whether statistics should be plotted by
  chromosome (see Details).

- stat.type:

  Character string of any combination of `"xo"` or `"dxo"` or `"miss"`.
  `"xo"` calculates the number of crossovers, `"dxo"` calculates the
  number of double crossover and `"miss"` calculates the number of
  missing values.

- id:

  Character string determining the column of `cross$pheno` that contains
  the genotype names.

- xo.lambda:

  A numerical value for the expected rate of recombination. (see
  Details).

- ...:

  Other arguments to be passed to the high level lattice plot.

## Details

This function uses `statGen` to profile statistics for the genotypes for
the current order of the linkage map. Any combination of `"xo"` or
`"dxo"` or `"miss"` may be given to simultaneous plot. If `bychr = TRUE`
then the plots will be further partitioned by linkage groups given by
`chr`.

If a numerical value is given for `xo.lambda` then the recombination
count for each genotype is tested against the expected recombination
rate `xo.lambda` using a simple one-tailed test of a Poisson mean. Any
lines that have a p-value less than than a family wise error rate based
on bonferroni adjustment of the usual alpha level of 0.05 are annotated
on the profiles being plotted.

## Value

A lattice panel plot with panels described by the `stat.type` given in
the call and genotype statistics are returned invisibly. If `xo.lambda`
is not NULL then these statistics also include a logical vector named
`"xo.lambda"` that is returned from testing the individuals for inflated
recombination rates (see Details).

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`statGen`](https://drj001.github.io/ASMap/reference/statGen.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## profile all genotype crossover and double crossover statistics

profileGen(mapDH, bychr = FALSE, stat.type = c("xo","dxo"),
     xo.lambda = 25, layout = c(1,3))

```
