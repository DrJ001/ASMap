# Push markers into an established R/qtl linkage map.

Push unlinked markers or markers that were originally placed aside by
`pullCross` back into linkage groups of an established R/qtl linkage
map.

## Usage

``` r
pushCross(object, type = c("co.located","seg.distortion",
          "missing","unlinked"), unlinked.chr = NULL,
          pars = NULL, ...)
```

## Arguments

- object:

  An R/qtl `cross` object with class structure `"bc"`, `"dh"`,
  `"riself"`, `"bcsft"`. (see
  [`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  for more details.)

- type:

  A character string determining the type of markers to be pushed into
  the linkage map (see Details).

- unlinked.chr:

  A character string of linkage group names containing markers that
  require pushing into the remaining linkage groups of the object. This
  is only useful when `type="unlinked"`. Default is `NULL`.

- pars:

  A list of parameters that are used by `pushCross` to push markers a
  certain type back into the linkage group. The default NULL calls the
  parameter initialization function `pp.init` with defaults (see Details
  and Examples).

- ...:

  Currently ignored.

## Details

This function was written explicitly to complement `pullCross` by
"pushing" markers of certain types back into linkage groups of an
established linkage map.

Currently supported marker types are:

- `type = "co.located"`. Users can push co-located markers back into the
  linkage map that have been set aside in the cross object element
  `co.located`. To ensure this can be used at any stage of the linkage
  map construction process the function disregards the linkage group
  information provided in the table formed by using `pullCross`. Instead
  it uses the current positions of the markers in the reduced linkage
  map to determine where to push the co-located markers back to.

- `type = "seg.distortion"`. Users can push markers from the
  `"seg.distortion"` element of the object back into a linkage map using
  the thresholding mechanisms `seg.thresh` and `seg.ratio` called using
  `pars`. If `seg.thresh` is given then the markers are pushed back that
  have p-values that are GREATER than `seg.thresh`. If `pars` contains
  an element `seg.ratio` then markers are pushed back based on the ratio
  provided. The ratio must be in character format and of the type
  "AA:BB" for two allele populations and "AA:AB:BB" for three allele
  populations (see Examples for more details). Markers are pushed back
  if their allele proportions are LESS than the largest proportional
  ratio or GREATER than the smallest proportional ratio given in
  `seg.thresh`. If neither thresholding mechanisms are given then the
  default is to use `seg.thresh = 0.05`.

- `type = "missing"`. Users can push markers from the object element
  `"missing"` back into the linkage map using the thresholding parameter
  `miss.thresh` called using `pars`. Markers will be pushed back that
  have a proportion of missing values LESS than `miss.thresh`. If no
  value is given for this parameter it defaults to 0.1 or 10% missing
  values.

- `type = "unlinked"`. Users can push unlinked markers that reside in
  linkage groups of the established linkage map. If this type is chosen
  `unlinked.chr` must be a character string of linkage group names in
  the object.

For types `"seg.distortion"`, `"missing"` and `"unlinked"` a fast
clustering method is used to allocate markers to established linkage
groups. This is done very efficiently by reducing the constructed
linkage map to a skeleton set of markers before checking linkages. How
these linkages are formed can be tweaked by setting `max.rf` and
`min.lod` when calling `pars`. These currently default to
`max.rf = 0.25` and `min.lod = 3`.

Users should explicitly avoid the use of "UL" as part of a linkage group
name as this is used internally to name unlinked groups of markers if
required. It should also be noted that this function does not
re-construct the object after allocating markers to linkage groups. For
efficient linkage map reconstruction of an R/qtl object see
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md).

## Value

The cross object is returned with an identical class structure as the
inputted cross object with additional markers from the marker types
pushed into linkage groups of the established linkage map. If all
markers of an element type are pushed back then the element type is
removed from the object.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`pullCross`](https://drj001.github.io/ASMap/reference/pullCross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## pull co-located markers from map

mapDH.c <- pullCross(mapDH, type = "co.located")
mapDH.c$co.located$table
#>    bins chr    mark
#> 1     1 1B1 1B1.m.4
#> 2     1 1B1 1B1.m.5
#> 3     2  1D  1D.m.3
#> 4     2  1D  1D.m.4
#> 5     3  2B 2B.m.31
#> 6     3  2B 2B.m.32
#> 7     4 2D1 2D1.m.6
#> 8     4 2D1 2D1.m.7
#> 9     5  3A 3A.m.35
#> 10    5  3A 3A.m.36
#> 11    6  3B 3B.m.11
#> 12    6  3B 3B.m.12
#> 13    7  3D  3D.m.1
#> 14    7  3D  3D.m.2
#> 15    8  4A  4A.m.7
#> 16    8  4A  4A.m.8
#> 17    9  5A 5A.m.36
#> 18    9  5A 5A.m.37
#> 19   10  5B 5B.m.15
#> 20   10  5B 5B.m.16
#> 21   11  5D  5D.m.4
#> 22   11  5D  5D.m.5
#> 23   12  6D  6D.m.2
#> 24   12  6D  6D.m.3
#> 25   13  7A 7A.m.32
#> 26   13  7A 7A.m.33
#> 27   14  7B 7B.m.21
#> 28   14  7B 7B.m.22
#> 29   15  7D  7D.m.2
#> 30   15  7D  7D.m.3

## push co-located markers back into linkage map

mapDH.z <- pushCross(mapDH.c, type = "co.located")
pull.map(mapDH.z)
#> $`1A`
#>     1A.m.1     1A.m.2     1A.m.3     1A.m.4     1A.m.5     1A.m.6     1A.m.7 
#>   0.000000   1.372970   4.126799  19.254021  20.629050  21.087686  21.546415 
#>     1A.m.8     1A.m.9    1A.m.10    1A.m.11    1A.m.12    1A.m.13    1A.m.14 
#>  22.005144  22.463637  25.219148  26.136448  27.513085  27.971731  27.971731 
#>    1A.m.15    1A.m.16    1A.m.17    1A.m.18    1A.m.19    1A.m.20    1A.m.21 
#>  28.430376  29.807059  30.265649  31.183296  31.641905  33.018542  33.936097 
#>    1A.m.22    1A.m.23    1A.m.24    1A.m.25    1A.m.26    1A.m.27    1A.m.28 
#>  34.394799  34.853244  38.069164  39.442113  68.470093  70.302626  70.760930 
#>    1A.m.29    1A.m.30    1A.m.31    1A.m.32    1A.m.33    1A.m.34    1A.m.35 
#>  73.976804  75.812508  76.122866  76.433365  76.892109  77.350044  85.190751 
#>    1A.m.36    1A.m.37    1A.m.38    1A.m.39    1A.m.40    1A.m.41 
#>  89.768584  91.144156  98.069974 100.364864 100.823405 101.272222 
#> 
#> $`1B1`
#>   1B1.m.1   1B1.m.2   1B1.m.3   1B1.m.4   1B1.m.5 
#> 0.0000000 0.9075244 3.6626600 7.7914418 7.7914418 
#> 
#> $`1B2`
#>    1B2.m.1    1B2.m.2    1B2.m.3    1B2.m.4    1B2.m.5    1B2.m.6    1B2.m.7 
#>  0.0000000  0.9069125  9.2418212  9.6996066 11.0761971 12.4527876 12.9112832 
#>    1B2.m.8    1B2.m.9   1B2.m.10   1B2.m.11   1B2.m.12   1B2.m.13   1B2.m.14 
#> 14.7471618 16.1237055 16.5823413 17.0410435 17.9586450 18.8762388 18.8762389 
#>   1B2.m.15   1B2.m.16   1B2.m.17   1B2.m.18   1B2.m.19   1B2.m.20   1B2.m.21 
#> 19.3350264 19.7936622 21.1703450 21.6286474 25.3055230 27.1410256 28.9760369 
#>   1B2.m.22   1B2.m.23   1B2.m.24   1B2.m.25   1B2.m.26   1B2.m.27   1B2.m.28 
#> 37.1409287 38.4752454 39.3926924 40.7693293 41.2279261 42.0729242 42.9177826 
#>   1B2.m.29   1B2.m.30   1B2.m.31   1B2.m.32   1B2.m.33 
#> 44.7537146 45.2122569 46.1297924 47.5054672 56.8652986 
#> 
#> $`1D`
#>    1D.m.1    1D.m.2    1D.m.3    1D.m.4    1D.m.5    1D.m.6    1D.m.7    1D.m.8 
#>  0.000000  2.285816  2.742698  2.742698 21.465574 25.652893 57.200980 58.104090 
#>    1D.m.9   1D.m.10 
#> 58.568087 81.061896 
#> 
#> $`2A`
#>    2A.m.1    2A.m.2    2A.m.3    2A.m.4    2A.m.5    2A.m.6    2A.m.7    2A.m.8 
#>   0.00000  14.60776  24.35533  25.27181  27.10621  41.22918  46.75093  47.20874 
#>    2A.m.9   2A.m.10   2A.m.11   2A.m.12   2A.m.13   2A.m.14   2A.m.15   2A.m.16 
#>  51.34564  67.48319  67.94012  68.39885  68.85758  69.31635  69.37518  71.31173 
#>   2A.m.17   2A.m.18   2A.m.19   2A.m.20   2A.m.21   2A.m.22   2A.m.23   2A.m.24 
#>  73.67100  75.50543  88.15661  88.61360  92.75224  93.21058  93.66903  96.88485 
#>   2A.m.25   2A.m.26   2A.m.27   2A.m.28   2A.m.29   2A.m.30   2A.m.31   2A.m.32 
#>  99.17995 101.52377 102.94898 103.86646 104.80807 120.36287 120.81940 125.88349 
#>   2A.m.33   2A.m.34   2A.m.35   2A.m.36   2A.m.37   2A.m.38   2A.m.39   2A.m.40 
#> 126.34122 131.71419 143.11651 144.22823 144.88222 146.25907 146.99787 149.16025 
#> 
#> $`2B`
#>    2B.m.1    2B.m.2    2B.m.3    2B.m.4    2B.m.5    2B.m.6    2B.m.7    2B.m.8 
#>  0.000000  2.776311 25.842944 26.298448 31.825602 41.104786 49.437816 51.272946 
#>    2B.m.9   2B.m.10   2B.m.11   2B.m.12   2B.m.13   2B.m.14   2B.m.15   2B.m.16 
#> 51.731535 52.189561 60.997121 61.913709 63.749541 65.126085 65.584721 66.043450 
#>   2B.m.17   2B.m.18   2B.m.19   2B.m.20   2B.m.21   2B.m.22   2B.m.23   2B.m.24 
#> 66.502132 67.419714 68.337204 69.713841 70.172383 71.548974 72.925564 73.384060 
#>   2B.m.25   2B.m.26   2B.m.27   2B.m.28   2B.m.29   2B.m.30   2B.m.31   2B.m.32 
#> 75.219928 75.219928 75.678650 76.596278 77.054961 77.513690 77.972325 77.972325 
#>   2B.m.33   2B.m.34   2B.m.35 
#> 79.350096 79.350096 95.477694 
#> 
#> $`2D1`
#>    2D1.m.1    2D1.m.2    2D1.m.3    2D1.m.4    2D1.m.5    2D1.m.6    2D1.m.7 
#>  0.0000000  0.4488416  0.9075944  1.3662863  2.3836941  5.9508951  5.9508951 
#>    2D1.m.8 
#> 19.1368202 
#> 
#> $`2D2`
#>    2D2.m.1    2D2.m.2    2D2.m.3    2D2.m.4    2D2.m.5    2D2.m.6    2D2.m.7 
#>  0.0000000  0.4488151  0.9317667  8.7708039 10.9066024 13.0424288 13.5006912 
#>    2D2.m.8    2D2.m.9   2D2.m.10   2D2.m.11   2D2.m.12   2D2.m.13 
#> 17.6783731 49.1819741 49.6150457 49.9659711 56.7376665 59.9433335 
#> 
#> $`3A`
#>      3A.m.1      3A.m.2      3A.m.3      3A.m.4      3A.m.5      3A.m.6 
#>   0.0000000   0.4488176   0.9068564   7.8339385  11.5115292  12.4281903 
#>      3A.m.7      3A.m.8      3A.m.9     3A.m.10     3A.m.11     3A.m.12 
#>  18.8870956  19.8017237  54.1887247  55.3015873  55.9562210  58.0853345 
#>     3A.m.13     3A.m.14     3A.m.15     3A.m.16     3A.m.17     3A.m.18 
#>  62.5185226  62.9768368  63.4355193  64.3531474  64.8121744  64.8121745 
#>     3A.m.19     3A.m.20     3A.m.21     3A.m.22     3A.m.23     3A.m.24 
#>  79.4373752  82.6518066  83.2246888  83.7972184  90.2562018  90.7137076 
#>     3A.m.25     3A.m.26     3A.m.27     3A.m.28     3A.m.29     3A.m.30 
#>  96.7061154 103.1534712 109.1380198 120.8185483 131.5257233 139.7732679 
#>     3A.m.31     3A.m.32     3A.m.33     3A.m.34     3A.m.35     3A.m.36 
#> 141.6944239 142.6117164 144.9072399 145.3657173 146.7424296 146.7424296 
#>     3A.m.37 
#> 147.1911542 
#> 
#> $`3B`
#>      3B.m.1      3B.m.2      3B.m.3      3B.m.4      3B.m.5      3B.m.6 
#>   0.0000000   0.4487713   1.3663994   1.8250820   2.2836236   4.5703638 
#>      3B.m.7      3B.m.8      3B.m.9     3B.m.10     3B.m.11     3B.m.12 
#>  29.2532005  35.8080316  36.7227763  58.6806410  59.1368252  59.1368252 
#>     3B.m.13     3B.m.14     3B.m.15     3B.m.16     3B.m.17     3B.m.18 
#>  59.5955505  60.0538965  64.1925357  64.6508816  65.1096201  65.1096202 
#>     3B.m.19     3B.m.20     3B.m.21     3B.m.22     3B.m.23     3B.m.24 
#>  65.5683123  66.4854592  71.5490426  76.5820437  79.3776517  93.9256998 
#>     3B.m.25     3B.m.26     3B.m.27     3B.m.28     3B.m.29     3B.m.30 
#>  98.0628133  98.2732301 102.1480203 109.9968916 122.1555489 132.3762327 
#> 
#> $`3D`
#>    3D.m.1    3D.m.2    3D.m.3    3D.m.4    3D.m.5    3D.m.6 
#> 0.0000000 0.0000000 0.4531971 3.6430045 8.2177422 8.6663099 
#> 
#> $`4A`
#>     4A.m.1     4A.m.2     4A.m.3     4A.m.4     4A.m.5     4A.m.6     4A.m.7 
#>   0.000000   1.332938   3.741521   3.888260   4.347020   4.805749   5.264338 
#>     4A.m.8     4A.m.9    4A.m.10    4A.m.11    4A.m.12    4A.m.13    4A.m.14 
#>   5.264338   7.091815  49.296209  51.581339  52.498880  52.956766  60.820566 
#>    4A.m.15    4A.m.16    4A.m.17    4A.m.18    4A.m.19    4A.m.20    4A.m.21 
#>  63.115352  63.570243  93.879162  94.334155  95.632810  99.231442 100.148746 
#>    4A.m.22    4A.m.23    4A.m.24    4A.m.25    4A.m.26    4A.m.27    4A.m.28 
#> 100.607383 101.524965 102.442547 102.901138 104.259240 120.373218 128.698538 
#>    4A.m.29    4A.m.30 
#> 140.366059 147.276666 
#> 
#> $`4B`
#>    4B.m.1    4B.m.2    4B.m.3    4B.m.4    4B.m.5    4B.m.6    4B.m.7    4B.m.8 
#>   0.00000  17.75680  18.21346  25.60830  28.36290  31.10695  45.22634  47.52211 
#>    4B.m.9   4B.m.10   4B.m.11   4B.m.12   4B.m.13   4B.m.14   4B.m.15   4B.m.16 
#>  47.98062  48.89827  49.35695  49.81568  50.27441  50.73314  51.19173  53.02769 
#>   4B.m.17   4B.m.18   4B.m.19   4B.m.20   4B.m.21   4B.m.22   4B.m.23   4B.m.24 
#>  53.48609  55.78160  56.69904  57.15763  58.53427  59.45112  66.37712  67.29406 
#>   4B.m.25   4B.m.26   4B.m.27   4B.m.28   4B.m.29   4B.m.30   4B.m.31   4B.m.32 
#>  67.75236  71.87789  79.72830  85.25552  86.63163  87.51716 101.11891 109.48525 
#> 
#> $`4D`
#>     4D.m.1     4D.m.2     4D.m.3     4D.m.4     4D.m.5     4D.m.6 
#>  0.0000000  0.9119796  2.2948855  4.1307497  8.2701607 12.8620983 
#> 
#> $`5A`
#>     5A.m.1     5A.m.2     5A.m.3     5A.m.4     5A.m.5     5A.m.6     5A.m.7 
#>   0.000000   6.449477   6.907428   8.743399   9.201988   9.660670  10.577817 
#>     5A.m.8     5A.m.9    5A.m.10    5A.m.11    5A.m.12    5A.m.13    5A.m.14 
#>  15.641291  21.699011  51.219400  51.650427  53.485972  58.086437  62.224519 
#>    5A.m.15    5A.m.16    5A.m.17    5A.m.18    5A.m.19    5A.m.20    5A.m.21 
#>  64.059271  72.394000  73.769747  75.146367  75.605032  76.063278  81.125938 
#>    5A.m.22    5A.m.23    5A.m.24    5A.m.25    5A.m.26    5A.m.27    5A.m.28 
#>  94.759321  95.674890 101.667865 102.125675 105.802551 107.638129 107.638130 
#>    5A.m.29    5A.m.30    5A.m.31    5A.m.32    5A.m.33    5A.m.34    5A.m.35 
#> 108.555730 110.391734 110.850277 111.767670 114.522902 117.738637 118.655974 
#>    5A.m.36    5A.m.37    5A.m.38    5A.m.39    5A.m.40    5A.m.41    5A.m.42 
#> 118.655974 118.655974 119.573455 121.409422 121.672425 126.538352 127.914577 
#>    5A.m.43    5A.m.44    5A.m.45    5A.m.46    5A.m.47    5A.m.48    5A.m.49 
#> 128.373213 128.831942 129.290671 129.748495 138.555758 140.390835 140.847232 
#>    5A.m.50    5A.m.51    5A.m.52    5A.m.53    5A.m.54 
#> 160.094028 177.720807 178.636839 179.095521 179.544339 
#> 
#> $`5B`
#>      5B.m.1      5B.m.2      5B.m.3      5B.m.4      5B.m.5      5B.m.6 
#>   0.0000000   0.4506512   2.7568703   3.6744561   4.5918528   6.8687311 
#>      5B.m.7      5B.m.8      5B.m.9     5B.m.10     5B.m.11     5B.m.12 
#>  16.6069342  17.0642166  21.6647206  26.7281663  29.4832023  30.3991798 
#>     5B.m.13     5B.m.14     5B.m.15     5B.m.16     5B.m.17     5B.m.18 
#>  43.5404968  43.9961988  58.6230583  58.6230583  59.0798487  62.7568481 
#>     5B.m.19     5B.m.20     5B.m.21     5B.m.22     5B.m.23     5B.m.24 
#>  62.7568482  63.6747200  65.0511340  68.2665322  73.7942486  75.1704016 
#>     5B.m.25     5B.m.26     5B.m.27     5B.m.28     5B.m.29     5B.m.30 
#>  75.6290375  76.0873834  80.2234080 103.2977100 103.7537667 103.7537668 
#>     5B.m.31     5B.m.32     5B.m.33     5B.m.34     5B.m.35     5B.m.36 
#> 104.6713214 105.5889034 106.0475396 106.9651677 107.4238503 107.8823913 
#>     5B.m.37     5B.m.38     5B.m.39     5B.m.40     5B.m.41     5B.m.42 
#> 110.1780101 110.6365936 111.5542804 112.0129826 112.4717114 112.9304403 
#>     5B.m.43     5B.m.44     5B.m.45     5B.m.46     5B.m.47     5B.m.48 
#> 113.3891425 114.3067902 114.7654727 115.2242016 115.6829305 116.1416594 
#>     5B.m.49     5B.m.50     5B.m.51     5B.m.52     5B.m.53     5B.m.54 
#> 116.6003883 117.0589769 118.8949474 119.3533482 121.6485891 125.3240331 
#>     5B.m.55     5B.m.56 
#> 138.4650325 139.3713760 
#> 
#> $`5D`
#>     5D.m.1     5D.m.2     5D.m.3     5D.m.4     5D.m.5     5D.m.6 
#>  0.0000000  0.4509275  0.9077973 17.5556736 17.5556736 18.0026320 
#> 
#> $`6A`
#>     6A.m.1     6A.m.2     6A.m.3     6A.m.4     6A.m.5     6A.m.6     6A.m.7 
#>  0.0000000  0.4508862  1.3683270  3.6637125  5.4993607  7.3158176 18.9666244 
#>     6A.m.8     6A.m.9    6A.m.10    6A.m.11    6A.m.12    6A.m.13    6A.m.14 
#> 24.0229475 24.4803012 33.1811682 34.4546825 42.2927177 46.4051258 53.3307152 
#>    6A.m.15    6A.m.16    6A.m.17    6A.m.18    6A.m.19    6A.m.20    6A.m.21 
#> 55.1659585 56.0834469 56.5421295 57.0008977 57.0008978 58.8366741 60.2124275 
#>    6A.m.22    6A.m.23    6A.m.24    6A.m.25    6A.m.26    6A.m.27 
#> 68.0701105 92.8558487 93.7705717 94.2290183 96.9845748 97.4331564 
#> 
#> $`6B`
#>      6B.m.1      6B.m.2      6B.m.3      6B.m.4      6B.m.5      6B.m.6 
#>   0.0000000   0.4509718   0.4509719   2.2867484   3.6632921   4.1209679 
#>      6B.m.7      6B.m.8      6B.m.9     6B.m.10     6B.m.11     6B.m.12 
#>  13.4002160  33.1827514  33.6392192  34.0974157  41.0089180  45.1327419 
#>     6B.m.13     6B.m.14     6B.m.15     6B.m.16     6B.m.17     6B.m.18 
#>  46.0499903  46.5085798  47.8852626  48.3438982  48.3438983  48.3438984 
#>     6B.m.19     6B.m.20     6B.m.21     6B.m.22     6B.m.23     6B.m.24 
#>  48.3438985  48.3438986  49.2614722  49.2614723  50.1789781  52.9346078 
#>     6B.m.25     6B.m.26     6B.m.27     6B.m.28     6B.m.29     6B.m.30 
#>  53.3923072  61.2308253  61.7139878  62.1727244  62.6312174  65.3848871 
#>     6B.m.31     6B.m.32     6B.m.33     6B.m.34     6B.m.35     6B.m.36 
#>  82.5455879  85.7602715  88.0555455  95.4495538  97.7443918  98.2028440 
#>     6B.m.37     6B.m.38     6B.m.39     6B.m.40     6B.m.41 
#>  99.1012419 100.9283497 103.5005078 104.6956501 107.9022594 
#> 
#> $`6D`
#>     6D.m.1     6D.m.2     6D.m.3     6D.m.4     6D.m.5     6D.m.6     6D.m.7 
#>  0.0000000  0.4509325  0.4509325  0.9096708  0.9096709  1.3684093  1.8260662 
#>     6D.m.8     6D.m.9    6D.m.10    6D.m.11    6D.m.12    6D.m.13    6D.m.14 
#> 12.0607962 12.5186442 13.4341937 31.8613959 60.6944586 67.1358830 70.3317538 
#>    6D.m.15 
#> 71.2434628 
#> 
#> $`7A`
#>     7A.m.1     7A.m.2     7A.m.3     7A.m.4     7A.m.5     7A.m.6     7A.m.7 
#>   0.000000   0.455080   2.299092   6.899104  15.233181  20.293404  46.857304 
#>     7A.m.8     7A.m.9    7A.m.10    7A.m.11    7A.m.12    7A.m.13    7A.m.14 
#>  47.312644  50.068210  50.418629  52.146105  53.522606  54.520517  59.039901 
#>    7A.m.15    7A.m.16    7A.m.17    7A.m.18    7A.m.19    7A.m.20    7A.m.21 
#>  60.416164  60.874753  61.792354  61.792354  62.251056  63.166222  84.570224 
#>    7A.m.22    7A.m.23    7A.m.24    7A.m.25    7A.m.26    7A.m.27    7A.m.28 
#>  94.324162  95.241500  97.535316 113.165545 121.027822 121.944226 126.545563 
#>    7A.m.29    7A.m.30    7A.m.31    7A.m.32    7A.m.33 
#> 127.004219 127.921563 131.152392 131.603038 131.603038 
#> 
#> $`7B`
#>     7B.m.1     7B.m.2     7B.m.3     7B.m.4     7B.m.5     7B.m.6     7B.m.7 
#>   0.000000   9.274417  15.269479  15.727100  21.137582  23.030407  23.200289 
#>     7B.m.8     7B.m.9    7B.m.10    7B.m.11    7B.m.12    7B.m.13    7B.m.14 
#>  23.889280  24.806885  25.265568  25.724250  26.641832  27.559227  29.854751 
#>    7B.m.15    7B.m.16    7B.m.17    7B.m.18    7B.m.19    7B.m.20    7B.m.21 
#>  30.313245  31.230873  31.689595  31.689595  33.525231  36.280604  37.197714 
#>    7B.m.22    7B.m.23    7B.m.24    7B.m.25    7B.m.26    7B.m.27    7B.m.28 
#>  37.197714  40.413680  40.871888  43.627166  46.842903  47.736943  66.432981 
#>    7B.m.29    7B.m.30    7B.m.31    7B.m.32    7B.m.33    7B.m.34    7B.m.35 
#>  67.348263  70.103217  75.931239  76.228114  76.686766  78.063118  81.740125 
#>    7B.m.36    7B.m.37 
#>  82.195836 105.826638 
#> 
#> $`7D`
#>   7D.m.1   7D.m.2   7D.m.3   7D.m.4 
#>  0.00000 16.19294 16.19294 21.69733 
#> 
```
