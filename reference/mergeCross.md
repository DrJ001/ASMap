# Merge linkage groups of an qtl cross object

Merges linkage groups of an qtl cross object from a user specified list.

## Usage

``` r
mergeCross(cross, merge = NULL, gap = 5)
```

## Arguments

- cross:

  An qtl `cross` object with any class structure.

- merge:

  A list with elements containing the linkage groups to be merged with
  each element named by the proposed linkage group name (see Examples).

- gap:

  The cM gap to put between the merged map elements in the complete
  linkage group.

## Details

This merging function allows you to perform multiple merges of two or
more linkage groups in one call. Users should ensure linkage group names
are correct and that proposed linkage group names do not already exist.

## Value

The cross object is returned with identical class structure as the
inputted cross object. The `"geno"` element should now contain merged
linkage groups for the user defined merges.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`breakCross`](https://drj001.github.io/ASMap/reference/breakCross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

mapDH1 <- breakCross(mapDH, split = list("4A" = "4A.m.8"))
pull.map(mapDH1)[["4A.1"]]
#>   4A.m.1   4A.m.2   4A.m.3   4A.m.4   4A.m.5   4A.m.6   4A.m.7   4A.m.8 
#> 0.000000 1.332938 3.741521 3.888260 4.347020 4.805749 5.264338 5.264338 
#> attr(,"class")
#> [1] "A"
pull.map(mapDH1)[["4A.2"]]
#>    4A.m.9   4A.m.10   4A.m.11   4A.m.12   4A.m.13   4A.m.14   4A.m.15   4A.m.16 
#>   0.00000  42.20439  44.48952  45.40706  45.86495  53.72875  56.02354  56.47843 
#>   4A.m.17   4A.m.18   4A.m.19   4A.m.20   4A.m.21   4A.m.22   4A.m.23   4A.m.24 
#>  86.78735  87.24234  88.54099  92.13963  93.05693  93.51557  94.43315  95.35073 
#>   4A.m.25   4A.m.26   4A.m.27   4A.m.28   4A.m.29   4A.m.30 
#>  95.80932  97.16742 113.28140 121.60672 133.27424 140.18485 
#> attr(,"class")
#> [1] "A"

mapDH2 <- mergeCross(mapDH1, merge = list("4A" = c("4A.1","4A.2")))
pull.map(mapDH2)[["4A"]]
#>     4A.m.1     4A.m.2     4A.m.3     4A.m.4     4A.m.5     4A.m.6     4A.m.7 
#>   0.000000   1.332938   3.741521   3.888260   4.347020   4.805749   5.264338 
#>     4A.m.8     4A.m.9    4A.m.10    4A.m.11    4A.m.12    4A.m.13    4A.m.14 
#>   5.264338  10.264338  52.468731  54.753861  55.671402  56.129289  63.993089 
#>    4A.m.15    4A.m.16    4A.m.17    4A.m.18    4A.m.19    4A.m.20    4A.m.21 
#>  66.287874  66.742765  97.051684  97.506678  98.805332 102.403964 103.321269 
#>    4A.m.22    4A.m.23    4A.m.24    4A.m.25    4A.m.26    4A.m.27    4A.m.28 
#> 103.779905 104.697487 105.615069 106.073660 107.431763 123.545740 131.871060 
#>    4A.m.29    4A.m.30 
#> 143.538582 150.449188 
#> attr(,"class")
#> [1] "A"
```
