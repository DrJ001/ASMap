# An unconstructed marker set for a backcross barley population

An unconstructed marker set for a backcross barley population in the
form of an qtl object.

## Usage

``` r
data(mapBCu)
```

## Format

This data relates to an unconstructed version of
[`mapBC`](https://drj001.github.io/ASMap/reference/mapBC.md) and
consists of 3023 markers genotyped on 326 individuals with markers
randomly assorted on one large linkage group. The data is in qtl format
with a class structure `c("bc","cross")`. See
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html)
documentation for more details on the format of this object. This data
set forms the basis of the worked example in Chapter 3 of the vignette
(see vignette("ASMap") for complete details)

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Examples

``` r

data(mapBCu, package = "ASMap")
```
