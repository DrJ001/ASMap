# Profile individual marker and interval statistics for an R/qtl cross object

Graphically profile individual marker and interval statistics for an
R/qtl cross object

## Usage

``` r
profileMark(cross, chr, stat.type = "marker", use.dist = TRUE,
          map.function = "kosambi", crit.val = NULL,
          display.markers = FALSE, mark.line = FALSE, ...)
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

- stat.type:

  Character string of either `"marker"` or `"interval"` or both. Also
  this can be a set of character strings relating to individual marker
  or interval statistics that want to be viewed simultaneously (see
  Details).

- use.dist:

  Logical value determining whether the actual map distances should be
  use to represent marker positions. If `FALSE` then markers are placed
  equidistant from each other.

- map.function:

  Character string of either `"kosambi"`, `"haldane"`, `"morgan"` or
  `"cf"` defining the map function to be used for interval related
  statistics.

- crit.val:

  The critical value to be used in displaying marker or intervals above
  a certain threshold (see Details).

- display.markers:

  A logical value determining whether marker names should be displayed
  on the bottom axis.

- mark.line:

  A logical value determining whether vertical lines should be drawn at
  marker positions. This may be useful to line up marker positions
  across several plots.

- ...:

  Other arguments to be passed to the high level lattice plot.

## Details

This graphical function calls the function `statMark` to retrieve marker
and interval statistics. If `"marker"` is given as the `stat.type` then
the complete set of marker statistics is plotted simultaneously. If
`"interval"` is given as the `stat.type` then the function
simultaneously plots the complete set of interval statistics. Both can
also be chosen.

This function also allows users to choose any combination of marker or
interval statistics they would like to view. The set of available marker
statistics that can be profiled are given below

- `"seg.dist"`: Profile the -log10 p-value. results from a test of
  segregation distortion for each marker.

- `"miss"`: Profile the proportion of missing values for each marker.

- `"prop"`: Profile the allele proportions for each marker.

- `"dxo"`: Profile the number of double crossovers occurring at each
  marker.

The set of available interval statistics that can be profiled are given
below

- `"erf"`: Profile the recombination fractions for the intervals.

- `"lod"`: Profile the LOD score for the test of no linkage between
  markers in an interval.

- `"dist"`: Profile the interval map distance taken from the map
  component of each linkage group.

- `"mrf"`: Profile the map recombination fraction for the intervals.

- `"recomb"`: Profile the actual number of recombinations within each of
  the intervals.

If `crit.val="bonf"` and marker statistics are plotted then any markers
that have p-value for the test of segregation distortion less than the
family wise error rate based on a bonferroni adjustment of the usual
0.05 alpha level, are annotated on each of the marker plots. If any
interval statistics are being plotted then any intervals that have a
p-value for the test of no linkage that is less than a bonferroni
adjustment of the usual 0.05 alpha level are annotated on each of the
interval statistics plots.

## Value

A lattice panel plot is displayed with panels described by the
`stat.type` given in the call and the complete marker/interval
statistics are returned invisibly. If `crit.val` is not NULL then both
the marker/interval statistics are returned with an extra logical column
called `"crit.val"` from testing markers for segregation distortion and
intervals for weak linkage (see Details).

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

`profileMark`

## Examples

``` r

data(mapDH, package = "ASMap")

## profile chosen statistics

profileMark(mapDH, stat.type = c("seg.dist","prop","erf"), layout =
      c(1,4), type = "l")
#> Warning: Some markers at the same position on chr 1B1,1D,2B,2D1,3A,3B,3D,4A,5A,5B,5D,6D,7A,7B,7D; use jittermap().

```
