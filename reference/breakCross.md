# Break linkage groups of an qtl cross object

Breaks linkage groups of an qtl cross object from a user specified list.

## Usage

``` r
breakCross(cross, split = NULL, suffix = "numeric", sep = ".")
```

## Arguments

- cross:

  An qtl `cross` object with any class structure.

- split:

  A list named by the linkage groups required for splitting and
  containing marker names immediately preceding where the splits are to
  be made (see Details).

- suffix:

  This can be a vector of character strings containing `"numeric"` or
  `"alpha"` specifying whether integers or letters are to be appended to
  the old linkage group names to form new names. This argument may also
  be list with elements named by the linkage groups that are in `split`
  and containing the new names for the split linkage groups (see
  Examples).

- sep:

  The character separator to be used to separate the linkage group name
  and the suffix.

## Details

The splitting of any linkage group only needs to be defined by the
markers immediately preceding where the splits are to be made. Multiple
splits in the one linkage group are possible as well as splitting across
multiple linkage groups with one call.

## Value

The cross object is returned with identical class structure as the
inputted cross object. The `"geno"` element will contain separate
linkage groups for the user defined splits.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`mergeCross`](https://drj001.github.io/ASMap/reference/mergeCross.md)

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

## manually choose suffix

mapDH1 <- breakCross(mapDH, split = list("4A" = "4A.m.8"),
                     suffix = list("4A" = c("4AA","4AB")))
```
