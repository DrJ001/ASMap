# Individual genotype statistics for an R/qtl cross object

Individual genotype statistics for the current linkage map order of and
R/qtl cross object

## Usage

``` r
statGen(cross, chr, bychr = TRUE, stat.type = c("xo","dxo",
        "miss"), id = "Genotype")
```

## Arguments

- cross:

  An R/qtl `cross` object with class structure `"bc"`, `"dh"`,
  `"riself"`, `"bcsft"`. (see
  [`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  for more details.)

- chr:

  Character vector of linkage group names used for subsetting the
  linkage map.

- bychr:

  Logical vector determining whether statistics should be plotted by
  chromosome (see Details).

- stat.type:

  Character string of any combination of `"xo"` or `"dxo"` or both.
  `"miss"`. `"xo"` calculates the number of crossovers, `"dxo"`
  calculates the number of double crossover and `"miss"` calculates the
  number of missing values.

- id:

  Character string determining the column of `cross$pheno` that contains
  the genotype names.

## Details

This function is used in `profileGen` to plot any combination of
returned linkage map statistics on a single graphical display.

## Value

A list with elements named by the `stat.type` used in the call. If
`bychr = TRUE` then each element is a data frame of statistics with
columns named by the linkage groups. If `bychr = FALSE` then each
element is a vector of statistics named by the `stat.type`.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`profileGen`](https://drj001.github.io/ASMap/reference/profileGen.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## produce all genotype crossover and double crossover statistics

sg <- statGen(mapDH, stat.type = c("xo","dxo"))
```
