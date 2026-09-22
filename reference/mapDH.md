# A constructed linkage map for a doubled haploid wheat population

A constructed linkage map for a doubled haploid wheat population in the
form of a constructed qtl object.

## Usage

``` r
data(mapDH)
```

## Format

This data relates to a fully constructed linkage map of 599 markers
genotyped on 218 individuals. The linkage map consists of 23 linkage
groups spanning the whole genome. 584 markers are from the orignal map
with an additonal 12 co-located markers and 3 slightly distorted
markers. The map was constructed using the MSTmap algorithm integrated
in[`mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
with geentic distances estimated using the `"kosambi"` mapping function.
The data is in qtl format with a class structure `c("bc","cross")`. See
`read.cross` documentation for more details on the format of this
object.

## Examples

``` r

data(mapDH, package = "ASMap")
```
