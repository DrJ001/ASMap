# Pull markers from a linkage map.

Pull markers of a certain type from a linkage map and place them aside
in the R/qtl object and, if appropriate, keeping their connections with
the reduced linkage map.

## Usage

``` r
pullCross(object, chr, type = c("co.located","seg.distortion",
          "missing"), pars = NULL, replace = FALSE, ...)
```

## Arguments

- object:

  An qtl `cross` object with class structure `"bc"` `"dh"`, `"riself"`,
  `"bcsft"`. (see
  [`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  for more details.)

- chr:

  A character vector of linkage group names with which to subset the
  linkage map before pulling any markers.

- type:

  A character string determining the type of markers to be pulled from
  the map (see Details).

- pars:

  A list of parameters that are used by `pullCross` to pull markers of
  certain type. The default NULL calls the parameter initialization
  function `pp.init` with defaults (see Details and Examples).

- replace:

  A logical value determining whether the markers and summary of marker
  information that is pulled from the map replaces information that is
  already residing in the `type` element of the object.

- ...:

  Currently ignored.

## Details

This function gives users the ability to "pull" markers of several
different types from the linkage map and place them in appropriately
named elements of the cross object. These elements can be examined by
the user and can even be "pushed" back using the complementary command
`pushCross`.

Currently supported types are:

- `type = "co.located"`. This type gives the user the ability to reduce
  a linkage map to a unique set of markers for the purpose of efficient
  map construction. Co-located markers are pulled from the linkage map
  using the technology of `findDupMarkers` from the qtl package and
  places them aside in a separate list element called `"co.located"`.
  This element contains the removed marker data as well as a table that
  displays the connections between the co-located markers with markers
  that remain in the linkage map. If required, this table is used by
  `pushCross` to "push" the co-located markers back into the linkage
  map.

- `type = "seg.distortion"`. Users can pull markers with segregation
  distortion from a linkage map with two different thresholding
  mechanisms called using `pars`. If the list argument `pars` is used
  with an element called `seg.thresh` then markers are pulled from the
  map if the p-value from the test for segregation distortion is LESS
  than `seg.thresh`. Values of `seg.thresh` must be between 0 and 1. If
  `pars` contains an element `seg.ratio` then markers are pulled from
  the map based on the ratio provided. The ratio must be in character
  format and of the type "AA:BB" for two allele populations and
  "AA:AB:BB" for three allele populations (see Examples for more
  details). Markers are pulled if their allele proportions are GREATER
  than the largest proportional ratio or LESS than the smallest
  proportional ratio given in `seg.thresh`. If neither thresholding
  mechanisms are given then the default is to use `seg.thresh = 0.05`.
  If markers are found matching the above criteria they are pulled from
  the linkage map and placed aside in an element called
  `"seg.distortion"`. This element contains the removed distorted marker
  data as well as a table summarizing each of the markers. See examples
  below for more detail.

- `type = "missing"`. Users can pull markers with a proportional amount
  of missing allele scores. If `pars` contains an element `miss.thresh`
  then markers are pulled from the linkage map that have a proportion of
  missing values GREATER than `miss.thresh`. If no value is given for
  `miss.thresh` then it defaults to 0.1 or 10% missing values. If
  markers are found matching the above criteria they are pulled from the
  map and are placed aside in an separate list element called
  `"missing"`. This element contains the removed marker data as well as
  a table summarizing each of the markers. See examples below for more
  detail.

## Value

The cross object is returned with identical class structure as the
inputted cross object and an additional elements corresponding to the
marker types being pulled from the map.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`pushCross`](https://drj001.github.io/ASMap/reference/pushCross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## pull co-located markers from linkage map

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

## pull distorted markers from linkage map using seg.thresh

mapDH.s <- pullCross(mapDH, type = "seg.distortion", pars =
           list(seg.thresh = 0.05))
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

## pull distorted markers from linkage map using seg.ratio

mapDH.s <- pullCross(mapDH, type = "seg.distortion", pars =
            list(seg.ratio = "56:44"))
mapDH.s$seg.distortion$table
#>       mark chr         pos neglog10P     missing        AA        AB
#> 1  1A.m.26  1A  68.4700927  1.375015 0.000000000 0.5688073 0.4311927
#> 2  1A.m.27  1A  70.3026255  1.375015 0.000000000 0.5688073 0.4311927
#> 3  1A.m.28  1A  70.7609297  1.237267 0.000000000 0.5642202 0.4357798
#> 4  1A.m.29  1A  73.9768038  1.671789 0.000000000 0.5779817 0.4220183
#> 5  1A.m.30  1A  75.8125081  1.375015 0.000000000 0.5688073 0.4311927
#> 6  1A.m.32  1A  76.4333654  1.519835 0.000000000 0.5733945 0.4266055
#> 7  1A.m.33  1A  76.8921092  1.671789 0.000000000 0.5779817 0.4220183
#> 8  1A.m.34  1A  77.3500443  1.830929 0.000000000 0.5825688 0.4174312
#> 9  1A.m.35  1A  85.1907512  1.375015 0.000000000 0.5688073 0.4311927
#> 10 1A.m.36  1A  89.7685841  1.375015 0.000000000 0.5688073 0.4311927
#> 11 1A.m.37  1A  91.1441558  1.830929 0.000000000 0.5825688 0.4174312
#> 12 1A.m.38  1A  98.0699737  1.671789 0.000000000 0.5779817 0.4220183
#> 13 1A.m.39  1A 100.3648637  1.519835 0.000000000 0.5733945 0.4266055
#> 14 1A.m.40  1A 100.8234047  1.671789 0.000000000 0.5779817 0.4220183
#> 15 1A.m.41  1A 101.2722224  1.519835 0.000000000 0.5733945 0.4266055
#> 16  2A.m.1  2A   0.0000000  1.375015 0.000000000 0.5688073 0.4311927
#> 17  2A.m.9  2A  51.3456384  1.671789 0.000000000 0.5779817 0.4220183
#> 18 2A.m.20  2A  88.6135963  1.237267 0.000000000 0.5642202 0.4357798
#> 19 2A.m.26  2A 101.5237715  1.121967 0.018348624 0.5607477 0.4392523
#> 20 2A.m.29  2A 104.8080733  1.338322 0.032110092 0.5687204 0.4312796
#> 21 2D2.m.1 2D2   0.0000000  1.519835 0.000000000 0.5733945 0.4266055
#> 22 2D2.m.2 2D2   0.4488151  1.375015 0.000000000 0.5688073 0.4311927
#> 23 2D2.m.3 2D2   0.9317667  1.600692 0.004587156 0.5760369 0.4239631
#> 24 3A.m.16  3A  64.3531474  1.375015 0.000000000 0.5688073 0.4311927
#> 25 3A.m.17  3A  64.8121744  1.451708 0.004587156 0.5714286 0.4285714
#> 26 3A.m.18  3A  64.8121745  1.237267 0.000000000 0.5642202 0.4357798
#> 27 3A.m.19  3A  79.4373752  1.375015 0.000000000 0.5688073 0.4311927
#> 28 3A.m.20  3A  82.6518066  1.237267 0.000000000 0.5642202 0.4357798
#> 29 3A.m.22  3A  83.7972184  1.237267 0.000000000 0.5642202 0.4357798
#> 30 3B.m.15  3B  64.1925357  1.997307 0.000000000 0.5871560 0.4128440
#> 31 3B.m.16  3B  64.6508816  2.170970 0.000000000 0.5917431 0.4082569
#> 32 3B.m.17  3B  65.1096201  1.872327 0.027522936 0.5849057 0.4150943
#> 33 3B.m.18  3B  65.1096202  1.997307 0.000000000 0.5871560 0.4128440
#> 34 3B.m.19  3B  65.5683123  1.830929 0.000000000 0.5825688 0.4174312
#> 35 3B.m.20  3B  66.4854592  1.830929 0.000000000 0.5825688 0.4174312
#> 36 3B.m.23  3B  79.3776517  1.600692 0.004587156 0.5760369 0.4239631
#> 37 3B.m.25  3B  98.0628133  1.237267 0.000000000 0.5642202 0.4357798
#> 38  4B.m.1  4B   0.0000000  1.600692 0.004587156 0.5760369 0.4239631
#> 39  4B.m.2  4B  17.7567993  1.519835 0.000000000 0.5733945 0.4266055
#> 40  4B.m.3  4B  18.2134551  1.121967 0.018348624 0.5607477 0.4392523
#> 41  4B.m.5  4B  28.3629008  1.671789 0.000000000 0.5779817 0.4220183
#> 42  4B.m.6  4B  31.1069464  1.375015 0.000000000 0.5688073 0.4311927
#> 43  5A.m.3  5A   6.9074284  1.237267 0.000000000 0.4357798 0.5642202
#> 44 5A.m.25  5A 102.1256749  1.237267 0.000000000 0.4357798 0.5642202
#> 45 5B.m.31  5B 104.6713214  1.375015 0.000000000 0.5688073 0.4311927
#> 46 5B.m.32  5B 105.5889034  1.375015 0.000000000 0.5688073 0.4311927
#> 47 5B.m.33  5B 106.0475396  1.519835 0.000000000 0.5733945 0.4266055
#> 48 5B.m.34  5B 106.9651677  1.237267 0.000000000 0.5642202 0.4357798
#> 49 5B.m.35  5B 107.4238503  1.375015 0.000000000 0.5688073 0.4311927
#> 50 5B.m.36  5B 107.8823913  1.237267 0.000000000 0.5642202 0.4357798
#> 51 5B.m.37  5B 110.1780101  1.175090 0.004587156 0.5622120 0.4377880
#> 52 5B.m.41  5B 112.4717114  1.237267 0.000000000 0.5642202 0.4357798
#> 53 5B.m.53  5B 121.6485891  1.237267 0.000000000 0.5642202 0.4357798
#> 54 6B.m.12  6B  45.1327419  1.237267 0.000000000 0.5642202 0.4357798
#> 55 6D.m.12  6D  60.6944586  1.756872 0.004587156 0.4193548 0.5806452
#> 56  7B.m.6  7B  23.0304066  1.769873 0.013761468 0.4186047 0.5813953
#> 57  7B.m.7  7B  23.2002891  1.309862 0.004587156 0.4331797 0.5668203
#> 58 7B.m.19  7B  33.5252305  1.237267 0.000000000 0.4357798 0.5642202
```
