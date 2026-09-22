# A constructed linkage map for a backcross barley population

A constructed linkage map for a backcross barley population in the form
of a constructed qtl object.

## Usage

``` r
data(mapBC)
```

## Format

This data relates to a fully constructed linkage map of 3019 markers
genotyped on 300 individuals spanning the 7 linkage groups of the barley
genome. The map was constructed using the MSTmap algorithm integrated in
[`mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
with geentic distances estimated using the `"kosambi"` mapping function.
The data is in qtl format with a class structure `c("bc","cross")`. See
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html)
documentation for more details on the format of this object.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Examples

``` r

data(mapBC, package = "ASMap")
```
