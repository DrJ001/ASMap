# Simulated constructed linkage map for a self pollinated F2 barley population

Simulated constructed linkage map for a self pollinated F2 barley
population in the form of an qtl object.

## Usage

``` r
data(mapF2)
```

## Format

This data relates to a fully constructed linkage map of 700 simulated
markers genotyped on 250 individuals. The map consists of 7 linkage
groups, each contaning 100 markers spanning an approximate linkage group
length of 200cM. The map was constructed using
[`mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
from the ASMap package and map distances were estimated using the
`"kosambi"` mapping function. The data is in R/qtl format with a class
structure `c("bcsft","cross")`.

## Examples

``` r

data(mapF2, package = "ASMap")
```
