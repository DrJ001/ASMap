# Pulling and pushing markers

## Pulling and pushing markers

Often in linkage map construction, some pruning of the markers occurs
before initial construction. This may be the removal of markers with a
proportion of missing values higher than some desired threshold as well
as markers that are significantly distorted from their expected
Mendelian segregation patterns. Often the removal is permanent and the
possible importance of some of these markers is overlooked. A preferable
system would be to identify and place the problematic markers aside with
the intention of checking their usefulness at a later stage of the
construction process.

The R/ASMap package contains two functions that allow you to do this.
The first is the function
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
which provide users with a mechanism to “pull” markers of certain types
from a linkage map and place them aside. The complementary function
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
then allows users to “push” markers back into the linkage map at any
stage of the construction process. These can be seen as helper functions
for more efficient construction of linkage maps (see
[`?pullCross`](https://drj001.github.io/ASMap/reference/pullCross.md)
and
[`?pushCross`](https://drj001.github.io/ASMap/reference/pushCross.md)
for complete details).

``` r

pullCross(object, chr, type = c("co.located","seg.distortion","missing"),
               pars = NULL, replace = FALSE, ...)
pushCross(object, chr, type = c("co.located","seg.distortion","missing","unlinked"),
               unlinked.chr = NULL, pars = NULL, replace = FALSE, ...)
pp.init(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh = 0.1, max.rf =
             0.25, min.lod = 3)
```

In the current version of the package, the helper functions share three
types of markers that can be “pulled/pushed” from linkage maps. These
include markers that are co-located with other markers, markers that
have some defined segregation distortion and markers with a defined
proportion of missing values. If the argument `type` is
`"seg.distortion"` or `"missing"` then the initialization function
`pp.init` is used to determine the appropriate threshold parameter
setting (`seg.thresh`, `seg.ratio`, `miss.thresh`) that will be used to
pull/push markers from the linkage map. Users can set their own
parameters by setting the `pars` argument (see examples below). For each
of the different types,
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
will pull markers from the map and place them in separate elements of
the cross object. Within the elements, vital information is kept that
can be accessed by
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
to push the markers back at a later stage of linkage map construction.

The function
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
also contains another marker type called `"unlinked"` which, in
conjunction with the argument `unlinked.chr`, allows users to push
markers from an unlinked linkage group in the `geno` element of the
`object` into established linkage groups. This mechanism becomes vital,
for example, when pushing new markers into an established linkage map.
An example of this is presented in section Unknown linkage groups.

Again, the constructed linkage map `mapDH` will be used to showcase the
power of
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
and
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md).
Markers are pulled from the map that are co-located with other markers,
have significant segregation distortion with a p-value less than 0.02
and have a missing value proportion greater than 0.03.

``` r

mapDHs <- pullCross(mapDH, type = "co.located")
mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.02))
mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.03))
names(mapDHs)
```

    [1] "geno"           "pheno"          "co.located"     "seg.distortion"
    [5] "missing"       

``` r

names(mapDHs$co.located)
```

    [1] "table" "data" 

The cross object now contains three new elements named by the marker
types that are pulled from the map. In each of the elements there are
two elements, a table of information for the markers that are pulled and
the actual marker data in genotype by marker format (i.e. exactly the
same as the data contained in the linkage groups themselves).

``` r

mapDHs$seg.distortion$table
```

          mark chr      pos neglog10P     missing        AA        AB
    1  1A.m.34  1A 77.35004  1.830929 0.000000000 0.5825688 0.4174312
    2  1A.m.37  1A 91.14416  1.830929 0.000000000 0.5825688 0.4174312
    3  3B.m.15  3B 64.19254  1.997307 0.000000000 0.5871560 0.4128440
    4  3B.m.16  3B 64.65088  2.170970 0.000000000 0.5917431 0.4082569
    5  3B.m.17  3B 65.10962  1.872327 0.027522936 0.5849057 0.4150943
    6  3B.m.18  3B 65.10962  1.997307 0.000000000 0.5871560 0.4128440
    7  3B.m.19  3B 65.56831  1.830929 0.000000000 0.5825688 0.4174312
    8  3B.m.20  3B 66.48546  1.830929 0.000000000 0.5825688 0.4174312
    9  6D.m.12  6D 60.69446  1.756872 0.004587156 0.4193548 0.5806452
    10  7B.m.6  7B 23.03041  1.769873 0.013761468 0.4186047 0.5813953

The table element of each of the marker types `"missing"` and
`"seg.distortion"` consists of a summary of positional information as
well as information from
[`geno.table()`](https://rdrr.io/pkg/qtl/man/geno.table.html) in the
R/qtl package.

``` r

head(mapDHs$co.located$table)
```

      bins chr    mark
    1    1 1B1 1B1.m.4
    2    1 1B1 1B1.m.5
    3    2  1D  1D.m.3
    4    2  1D  1D.m.4
    5    3  2B 2B.m.31
    6    3  2B 2B.m.32

The table element for the marker type `"co.located"` contains
information on the markers that are co-located and the group or bin they
belong to. In each group the first marker is the reference marker that
remains in the linkage map and the remaining markers are pulled from the
map and placed aside.

Suppose now the map is undergone a construction or re-construction
process so that the linkage groups are artificially renamed. To do this
we will re-run MSTmap with `bychr = FALSE`.

``` r

mapDHs <- mstmap(mapDHs, bychr = FALSE, dist.fun = "kosambi", trace = TRUE, anchor = TRUE)
nmar(mapDHs)
```

     L.1 L.10 L.11 L.12 L.13 L.14 L.15 L.16 L.17 L.18 L.19  L.2 L.20 L.21 L.22 L.23 
      38   23    5    8   21   32    6   52   55    5   27    4   41   13   32   34 
    L.24  L.3  L.4  L.5  L.6  L.7  L.8  L.9 
       3   31    9   37   34    7   13   35 

The markers can now be pushed back into the linkage map using
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md).
The complete set of co.located markers are pushed back as well as
markers that have significant segregation distortion with p-values
greater than 0.001 and markers that have a missing value proportion less
than 0.05.

``` r

mapDHs <- pushCross(mapDHs, type = "co.located")
mapDHs <- pushCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.001))
mapDHs <- pushCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.05))
names(mapDHs)
```

    [1] "geno"  "pheno"

With the above parameter settings all markers from each marker type are
pushed back into the map and the marker type elements are removed from
the `object`.

``` r

pull.map(mapDHs)[[4]]
```

      4A.m.1   4A.m.2   4A.m.4   4A.m.3   4A.m.5   4A.m.6   4A.m.7   4A.m.8 
    0.000000 1.376494 3.671683 4.130412 5.047946 5.506674 5.965403 5.965403 
      4A.m.9 
    7.801089 
    attr(,"class")
    [1] "A"

``` r

pull.map(mapDHs)[[21]]
```

       2B.m.1    2B.m.2    2B.m.3    2B.m.4    2B.m.5    2B.m.6    2B.m.7    2B.m.8 
     0.000000  2.779284 25.847324 26.306053 31.833042 41.112443 49.445632 51.281319 
       2B.m.9   2B.m.10   2B.m.11   2B.m.12   2B.m.13   2B.m.14   2B.m.15   2B.m.16 
    51.740047 52.198776 61.004291 61.921825 63.757512 65.134006 65.592735 66.051463 
      2B.m.17   2B.m.18   2B.m.19   2B.m.20   2B.m.21   2B.m.22   2B.m.23   2B.m.24 
    66.510191 67.427726 68.345260 69.721754 70.180483 71.556977 72.933472 73.392200 
      2B.m.25   2B.m.26   2B.m.27   2B.m.28   2B.m.29   2B.m.30   2B.m.31   2B.m.32 
    75.227887 75.686615 76.604149 77.521684 77.980412 78.439141 78.897869 78.897869 
      2B.m.34   2B.m.33   2B.m.35 
    80.274363 80.733092 97.376797 
    attr(,"class")
    [1] "A"

For co-located markers the reference marker for each group is used as a
guide to place the set of co-locating markers back into the linkage map
adjacent to the reference marker. For example, the co-locating marker
`1D.m.4` on 1D appears adjacent to its reference marker `1D.m.3`.
Markers from the `"seg.distortion"` and `"missing"` elements are pushed
back to the end of the linkage groups ready for the map to be
re-constructed by chromosome. For example the distorted marker `6D.m.12`
is pushed to the end of 6D. A final run of MSTmap within each linkage
group will produce the desired map with all markers in their optimal
position.

``` r

mapDHs <- mstmap(mapDHs, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, anchor = TRUE, p.value = 2)
```
