# Graphical linkage group identity and alignment.

A graphical tool for identity and alignment of linkage groups in qtl
cross objects using reference maps.

## Usage

``` r
alignCross(object, chr, maps, ...)
```

## Arguments

- object:

  A qtl cross `object` with any class structure.

- chr:

  A character string of linkage group names or a logical vector equal to
  the length of the number of linkage groups (see
  [`subset.cross`](https://rdrr.io/pkg/qtl/man/subset.cross.html)).

- maps:

  A named list of qtl cross objects or `data.frame` objects containing
  markers that are present in `object`. The matching markers are used to
  identify and orient the `object` linkage groups (see Details).

- ...:

  Other arguments to be passed to the high level lattice plot.

## Details

If any list elements of `map` are qtl `"cross"` objects then marker
names, linkage group identity and genetic distance information are
extracted. List elements of `map` that are `data.frame` objects must
explicitly contain named columns `"marker"`, `"ref.chr"`, `"ref.dist"`
otherwise an error will be produced.

For each linkage group determined by `chr`, the contents of the listed
`maps` are checked for matching markers in `object`. For each `chr` and
reference map combination, a scatter plot of the `object` genetic
distances against the reference distances is displayed with reference
linkage group names as the plotting character. If a linkage group is in
correct orientation the overall slope of the scatter plot should be
positive. If a linkage group requires inverting then the overall slope
should be negative.

## Value

A lattice panel plot is displayed with panels labelled by a combination
of `chr` and the `maps` used as a reference. A data frame of these
results is also invisibly returned.

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

chrl <- sample(c(TRUE,FALSE), 23, replace = TRUE)
mapDH1 <- subset(mapDH, chr = chrl)
alignCross(mapDH, maps = list(DH = mapDH1), layout = c(3,5), col = 1:7)

```
