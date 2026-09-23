# Manipulating linkage maps

## Miscellaneous additional R/qtl functions

### Breaking and merging linkage groups

During the linkage map construction process there may be a requirement
to break or merge linkage groups. R/ASMap provides two functions to
achieve this.

``` r

breakCross(cross, split = NULL, suffix = "numeric", sep = ".")
mergeCross(cross, merge = NULL, gap = 5)
```

The
[`breakCross()`](https://drj001.github.io/ASMap/reference/breakCross.md)
function allows users to break linkage groups in a variety ways. The
`split` argument takes a list with elements named by the linkage group
names that require splitting and containing the markers that immediately
proceed where the splits are to be made. For example, a split of the
linkage group 3B and 6A linkage map `mapDH` after the seventh and
fifteenth marker respectively can be easily made using

``` r

mapDHb1 <- breakCross(mapDH, split = list("3B" = "3B.m.7","6A" = "6A.m.15"))
nmar(mapDHb1)
```

      1A  1B1  1B2   1D   2A   2B  2D1  2D2   3A 3B.1 3B.2   3D   4A   4B   4D   5A 
      41    5   33   10   40   35    8   13   37    7   23    6   30   32    6   54 
      5B   5D 6A.1 6A.2   6B   6D   7A   7B   7D 
      56    6   15   12   41   15   33   37    4 

The `split` argument is flexible and can handle multiple linkage groups
as well as multiple markers within linkage groups. The default use of
the `suffix` argument produces a numerical suffix attachment to the
original linkage groups being split with a separated by `sep`. Users can
also provide their own complete names for the new split linkage groups
by explicitly naming them in the `suffix` argument.

``` r

mapDHb2 <- breakCross(mapDH, split = list("3B" = "3B.m.7"), suffix = list("3B" = c("3B1","3B2")))
nmar(mapDHb2)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A 3B1 3B2  3D  4A  4B  4D  5A  5B  5D  6A  6B 
     41   5  33  10  40  35   8  13  37   7  23   6  30  32   6  54  56   6  27  41 
     6D  7A  7B  7D 
     15  33  37   4 

The
[`mergeCross()`](https://drj001.github.io/ASMap/reference/mergeCross.md)
function provides a method for merging linkage groups. Its argument
`merge` requires a list with elements named by the proposed linkage
group names required and containing the linkage groups to be merged. For
the linkage map `mapDHb1` containing split linkage groups 3B and 6A
created by a call to
[`breakCross()`](https://drj001.github.io/ASMap/reference/breakCross.md)
the call to
[`mergeCross()`](https://drj001.github.io/ASMap/reference/mergeCross.md)
would be

``` r

mapDHm <- mergeCross(mapDHb1, merge = list("3B" = c("3B.1","3B.2"),"6A" = c("6A.1","6A.2")))
nmar(mapDHm)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D 
     41   5  33  10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15 
     7A  7B  7D 
     33  37   4 

It should be noted that this function places an artificial genetic
distance gap between the merged linkage groups set by the `gap`
argument. Accurate distance estimation would require a separate map
estimation procedure after the merging has taken place.

### Rapid genetic distance estimation

The linkage map estimation function in R/qtl called
[`est.map()`](https://rdrr.io/pkg/qtl/man/est.map.html) can be invoked
individually or can be applied through
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html) when
setting the argument `estimate.map = TRUE`. The function applies the
multi-locus hidden Markov model technology of ([Lander and Green
1987](#ref-lg87)) to perform its calculations. Unfortunately this
technology is computationally cumbersome if there are many markers on a
linkage group and becomes more so if there are many missing allele calls
and genotyping errors present.

R/ASMap contains a small map estimation function that circumvents this
computational burden.

``` r

quickEst(object, chr, map.function = "kosambi", ...)
```

The function makes use of another function in R/qtl called
[`argmax.geno()`](https://rdrr.io/pkg/qtl/man/argmax.geno.html). This
function is also a multi-locus hidden Markov algorithm that uses the
observed markers present in a linkage group to impute pseudo-markers at
any chosen cM genetic distance. In this case, we only require a
reconstruction or imputation at the markers themselves. For the most
accurate imputation to occur there needs to be an estimate of genetic
distance already in place. To obtain an initial estimate of distance
[`est.rf()`](https://rdrr.io/pkg/qtl/man/est.rf.html) is called within
each linkage group defined by `chr` and recombination fractions are
converted to genetic distances based on `map.function`. Unlike
[`est.map()`](https://rdrr.io/pkg/qtl/man/est.map.html) the
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md)
function lives up to its namesake by providing the quickest genetic
distance calculations for large linkage maps.

``` r

map1 <- est.map(mapDH, map.function = "kosambi")
map1 <- subset(map1, chr = names(nmar(map1))[6:15])
map2 <- quickEst(mapDH, map.function = "kosambi")
map2 <- subset(map2, chr = names(nmar(map2))[6:15])
plotMap(map1, map2)
```

![Comparison of \`mapDH\` using \`est.map\` and
\`quickEst\`.](manipulating-maps_files/figure-html/quick2-1.png)

Comparison of `mapDH` using `est.map` and `quickEst`.

The linkage map `mapDH` was re-estimated using
[`est.map()`](https://rdrr.io/pkg/qtl/man/est.map.html) and
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md) and
a comparison of the resulting maps are given in the figure below. The
graphic indicates that there negligible changes in marker placement and
overall linkage group distances between the two linkage maps.

### Subsetting in R/ASMap

The functions
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
and
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
described in section Pulling and pushing markers are used to create and
manipulate extra list elements `"co.located"`, `"seg.distortion"` and
`"missing"` associated with different marker types. Each element
contains a data element consisting of a marker matrix equivalent in row
dimension to the marker elements of the linkage map they were pulled
from. Unfortunately, these list elements are not recognized by the
native R/qtl functions. If the R/qtl function `subset.cross()` is used
to subset the object to a reduced number of individuals then the data
component of each of these elements will not be subsetted accordingly.
In addition, the statistics in the table component of the elements
`"seg.distortion"` and `"missing"` will be incorrect for the newly
subsetted linkage map.

The
[`subsetCross()`](https://drj001.github.io/ASMap/reference/subsetCross.md)
function available in the R/ASMap package contains identical
functionality to `subset.cross` but also ensures the data components of
the extra list elements `"co.located"`, `"seg.distortion"` and
`"missing"` are subsetted to match the linkage map. In addition, for
elements `"seg.distortion"` and `"missing"` it also updates the table
components to reflect the newly subsetted map. This update ensures that
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
uses the most accurate information when deciding which markers to push
back into the linkage map.

Using the default `seg.thresh = 0.05` for the linkage map `mapDH`,
distorted markers are pulled from the map

``` r

mapDH.s <- pullCross(mapDH, type = "seg.distortion")
mapDH.s <- subsetCross(mapDH.s, ind = 3:218)
dim(mapDH.s$seg.distortion$data)[1]
```

    [1] 216

In this example the use of
[`subsetCross()`](https://drj001.github.io/ASMap/reference/subsetCross.md)
ensures that the data component of the `"seg.distortion"` element is the
same dimension as the map. The table element is also updated to ensure
the statistics are correct for the reduced subset of lines.

### Combining maps

Over the period of time that I have been involved in linkage map
construction there has been many occasions where I have required a
function that could merge R/qtl cross objects together in an intelligent
manner. For example, the merging of two linkage maps from the same
population that were independently built with markers from two different
platforms. This idea was the motivation behind the
[`combineMap()`](https://drj001.github.io/ASMap/reference/combineMap.md)
function in the R/ASMap package. The aim of the function was to merge
linkage maps based on map information, readying the combined linkage
groups for reconstruction through an efficient linkage map construction
process such as
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md).

``` r

combineMap(..., id = "Genotype", keep.all = TRUE)
```

The function takes an unlimited number of maps through the `...`
argument. The linkage maps must all have the same cross class structure
and contain the same genotype identifier `id`. At the current stage of
writing this vignette the function required unique maker names across
all linkage maps. This is expected to be relaxed at a later date so
linkage maps that share markers, such as nested association mapping
populations, can be merged effectively.

The merging of the maps happens intelligently using several components
of the map. Firstly the linkage maps are merged based on commonality
between the genotypes. If `keep.all = TRUE` the new combined linkage map
is “padded out” with missing values where genotypes are not shared. If
`keep.all= FALSE` the combined map is reduced to genotypes that are
shared among all linkage maps. Secondly, if linkage group names are
shared between maps then the markers from the shared linkage groups are
merged. To exemplify its use a duplicate of `mapDH` is made and 10
linkage group names have been altered, with the marker names inside each
of the linkage groups also altered to ensure they are unique.

``` r

mapDH1 <- mapDH
names(mapDH1$geno)[5:14] <- paste("L",1:10, sep = "")
mapDH1$geno <- lapply(mapDH1$geno, function(el){
  names(el$map) <- dimnames(el$data)[[2]] <- paste(names(el$map), "A", sep = "")
  el})
mapDHc <- combineMap(mapDH, mapDH1)
nmar(mapDHc)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D 
     82  10  66  20  40  35   8  13  37  30   6  30  32   6 108 112  12  54  82  30 
     7A  7B  7D  L1  L2  L3  L4  L5  L6  L7  L8  L9 L10 
     66  74   8  40  35   8  13  37  30   6  30  32   6 

The resulting combined map includes a combined marker set for the
linkage groups that shared the same name and distinct linkage groups for
unshared names. Again, note the function only merges the linkage maps
and does not reconstruct the final combined linkage map.

The advantages of this function may not be obvious at a first glance. If
an attempt is made to completely reconstruct the super set of linkage
maps, rather than combine them first, the identification of linkage
groups is lost. This function serves to preserve the important identity
of linkage groups. More examples of its use in common map construction
procedures will be explored in section chapter Post-construction linkage
map development of the next chapter.

Lander, E. S, and P Green. 1987. “Construction of Multilocus Genetic
Linkage Maps in Humans.” *Proceedings of the National Academy of
Science* 84: 2363–67.
