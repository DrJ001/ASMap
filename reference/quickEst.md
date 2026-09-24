# Very quick estimation of genetic map distances.

Very quick estimation of genetic map distances for a constructed R/qtl
object

## Usage

``` r
quickEst(object, chr, map.function = "kosambi", ...)
```

## Arguments

- object:

  An R/qtl `object` object with any class structure.

- chr:

  A character string of linkage group names that require (re)estimation
  of their genetic map distances.

- map.function:

  Character string of either `"kosambi"`, `"haldane"`, `"morgan"` or
  `"cf"` defining the mapping function to be used.

- ...:

  Other arguments passed to `argmax.geno`.

## Details

For linkage groups with large numbers of markers, the Hidden Markov
algorithm in `est.map` can be extremely slow. The computational burden
for this algorithm increases as the number of missing values and
genotyping errors increase. `quickEst` circumvents this by using the
Viterbi algorithm computationally implemented in `argmax.geno` of the
qtl package. Initial conservative estimates of the map distances are
calculated from inverting recombination fractions outputted from
`est.rf`. These are then passed to `argmax.geno` and imputation of
missing allele scores is performed along with re-estimation of map
distances.

## Value

The cross object is returned with identical class structure as the
inputted cross object.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`est.map`](https://rdrr.io/pkg/qtl/man/est.map.html)

## Examples

``` r

data(mapDH, package = "ASMap")

mapDH1 <- quickEst(mapDH, map.function = "kosambi")
```
