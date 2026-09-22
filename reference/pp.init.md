# Parameter initialization function

Parameter initialization function for `pushCross` and `pullCross`

## Usage

``` r
pp.init(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh = 0.1,
        max.rf = 0.25, min.lod = 3)
```

## Arguments

- seg.thresh:

  Numerical value between zero and one determining the p-value threshold
  for the test of marker segregation distortion.

- seg.ratio:

  A character string of the form "AA:BB" or "AA:AB:BB" describing the
  ratio of the alleles.

- miss.thresh:

  Numerical value between zero and one determining the proportion of
  missing values.

- max.rf:

  The maximum recombination fraction to consider when attempting to
  cluster pushed markers back into linkage groups.

- min.lod:

  The minimum LOD score to consider when attempting to cluster pushed
  markers back into linkage groups.

## Details

This parameter initialization function is used by the function
`pullCross` to pull markers from a linkage map and `pushCross` to push
markers back into a linkage map. How the arguments `seg.thresh`,
`seg.ratio` and `miss.thresh` are used depends on which function is
called. See `pushCross` and `pullCross` for more details.

## Value

Return user defined parameter values for each of the parameters.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`pushCross`](https://drj001.github.io/ASMap/reference/pushCross.md);
[`pullCross`](https://drj001.github.io/ASMap/reference/pullCross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## pull markers from a linkage map with a segregation distortion

pars <- pp.init(seg.thresh = 0.05)
mapDH.s <- pullCross(mapDH, type = "seg.distortion", pars = pars)
mapDH.s$seg.distortion$table
#>       mark chr         pos neglog10P     missing        AA        AB
#> 1  1A.m.26  1A  68.4700927  1.375015 0.000000000 0.5688073 0.4311927
#> 2  1A.m.27  1A  70.3026255  1.375015 0.000000000 0.5688073 0.4311927
#> 3  1A.m.29  1A  73.9768038  1.671789 0.000000000 0.5779817 0.4220183
#> 4  1A.m.30  1A  75.8125081  1.375015 0.000000000 0.5688073 0.4311927
#> 5  1A.m.32  1A  76.4333654  1.519835 0.000000000 0.5733945 0.4266055
#> 6  1A.m.33  1A  76.8921092  1.671789 0.000000000 0.5779817 0.4220183
#> 7  1A.m.34  1A  77.3500443  1.830929 0.000000000 0.5825688 0.4174312
#> 8  1A.m.35  1A  85.1907512  1.375015 0.000000000 0.5688073 0.4311927
#> 9  1A.m.36  1A  89.7685841  1.375015 0.000000000 0.5688073 0.4311927
#> 10 1A.m.37  1A  91.1441558  1.830929 0.000000000 0.5825688 0.4174312
#> 11 1A.m.38  1A  98.0699737  1.671789 0.000000000 0.5779817 0.4220183
#> 12 1A.m.39  1A 100.3648637  1.519835 0.000000000 0.5733945 0.4266055
#> 13 1A.m.40  1A 100.8234047  1.671789 0.000000000 0.5779817 0.4220183
#> 14 1A.m.41  1A 101.2722224  1.519835 0.000000000 0.5733945 0.4266055
#> 15  2A.m.1  2A   0.0000000  1.375015 0.000000000 0.5688073 0.4311927
#> 16  2A.m.9  2A  51.3456384  1.671789 0.000000000 0.5779817 0.4220183
#> 17 2A.m.29  2A 104.8080733  1.338322 0.032110092 0.5687204 0.4312796
#> 18 2D2.m.1 2D2   0.0000000  1.519835 0.000000000 0.5733945 0.4266055
#> 19 2D2.m.2 2D2   0.4488151  1.375015 0.000000000 0.5688073 0.4311927
#> 20 2D2.m.3 2D2   0.9317667  1.600692 0.004587156 0.5760369 0.4239631
#> 21 3A.m.16  3A  64.3531474  1.375015 0.000000000 0.5688073 0.4311927
#> 22 3A.m.17  3A  64.8121744  1.451708 0.004587156 0.5714286 0.4285714
#> 23 3A.m.19  3A  79.4373752  1.375015 0.000000000 0.5688073 0.4311927
#> 24 3B.m.15  3B  64.1925357  1.997307 0.000000000 0.5871560 0.4128440
#> 25 3B.m.16  3B  64.6508816  2.170970 0.000000000 0.5917431 0.4082569
#> 26 3B.m.17  3B  65.1096201  1.872327 0.027522936 0.5849057 0.4150943
#> 27 3B.m.18  3B  65.1096202  1.997307 0.000000000 0.5871560 0.4128440
#> 28 3B.m.19  3B  65.5683123  1.830929 0.000000000 0.5825688 0.4174312
#> 29 3B.m.20  3B  66.4854592  1.830929 0.000000000 0.5825688 0.4174312
#> 30 3B.m.23  3B  79.3776517  1.600692 0.004587156 0.5760369 0.4239631
#> 31  4B.m.1  4B   0.0000000  1.600692 0.004587156 0.5760369 0.4239631
#> 32  4B.m.2  4B  17.7567993  1.519835 0.000000000 0.5733945 0.4266055
#> 33  4B.m.5  4B  28.3629008  1.671789 0.000000000 0.5779817 0.4220183
#> 34  4B.m.6  4B  31.1069464  1.375015 0.000000000 0.5688073 0.4311927
#> 35 5B.m.31  5B 104.6713214  1.375015 0.000000000 0.5688073 0.4311927
#> 36 5B.m.32  5B 105.5889034  1.375015 0.000000000 0.5688073 0.4311927
#> 37 5B.m.33  5B 106.0475396  1.519835 0.000000000 0.5733945 0.4266055
#> 38 5B.m.35  5B 107.4238503  1.375015 0.000000000 0.5688073 0.4311927
#> 39 6D.m.12  6D  60.6944586  1.756872 0.004587156 0.4193548 0.5806452
#> 40  7B.m.6  7B  23.0304066  1.769873 0.013761468 0.4186047 0.5813953
#> 41  7B.m.7  7B  23.2002891  1.309862 0.004587156 0.4331797 0.5668203
```
