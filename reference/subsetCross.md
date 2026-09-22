# Subset an R/qtl object

Subset an R/qtl object by chromosome or by individuals for populations
used within the R/ASMap package.

## Usage

``` r
subsetCross(cross, chr, ind, ...)
```

## Arguments

- cross:

  A `"cross"` object generated from the R/qtl package. Specifically the
  object needs to inherit from one of the following classes `"bc"`,
  `"dh"`, `"riself"`, `"bcsft"` (see Details).

- chr:

  Optional vector specifying which chromosomes to keep or discard. This
  may be a logical, numeric, or character string vector. See
  [`?subset.cross`](https://rdrr.io/pkg/qtl/man/subset.cross.html).

- ind:

  Optional vector specifying which individuals to keep or discard. This
  may be a logical or numeric vector (see Details).

- ...:

  Kept for compatability with `subset.cross` and is ignored at this
  point.

## Details

This function is a replacement version of `subset.cross` that should be
used if the `cross` object contains any or all of the components
`"co.located"`, `"seg.distortion"` and `"missing"` created by a
`pullCross` call. For a given `ind`, the function calls `subset.cross`
to ensure that all elements created from calls to native R/qtl functions
are subsetted appropriately. In addition, the `"co.located"`,
`"seg.distortion"` and `"missing"` elements are also subsetted and if
components `"seg.distortion"` and `"missing"` exist, statistics in their
respective tables are recalculated.

It provides identical functionality to `subset.cross` with the exception
that `ind` can only be a logical or numeric vector.

## Value

The cross object is returned with the appropriate subsetting.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`subset.cross`](https://rdrr.io/pkg/qtl/man/subset.cross.html) and
`pullCross`

## Examples

``` r

data(mapDH, package = "ASMap")

mapDH.s <- pullCross(mapDH, type = "seg.distortion")
mapDH.s <- subsetCross(mapDH.s, ind = 3:218)
dim(mapDH.s$seg.distortion$data)
#> [1] 216  41
```
