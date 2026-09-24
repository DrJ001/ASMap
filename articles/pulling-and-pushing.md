# Pulling and pushing markers

Linkage map construction is commonly preceded by some pruning of the
marker set, typically the removal of markers exhibiting a high
proportion of missing values or significant departure from their
expected Mendelian segregation ratios. Such removal is usually
permanent, and the potential value of the discarded markers is thereby
overlooked. A marker that appears unsuitable when judged against
unordered data may prove entirely well behaved once a map has been
established.

A preferable approach is to identify the problematic markers and place
them aside, with the intention of reassessing their usefulness at a
later stage of the construction process. The functions
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
and
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
provide this facility. Rather than deleting a marker,
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
transfers it to a holding element of the cross object, retaining the
information required to restore it;
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
subsequently reintroduces it into an established map.

The constructed map `mapDH` is used below to illustrate both operations.

## Marker types

Three types of marker may be set aside.

| `type` | Markers affected | Rationale |
|----|----|----|
| `"co.located"` | Markers separated from another marker by zero distance | They contribute no ordering information, yet add to the computational cost of construction |
| `"seg.distortion"` | Markers whose segregation distortion p-value falls below `seg.thresh` | Distorted markers may not map to a unique location |
| `"missing"` | Markers whose proportion of missing values exceeds `miss.thresh` | Sparsely scored markers are positioned less reliably |

The relevant thresholds are supplied by
[`pp.init()`](https://drj001.github.io/ASMap/reference/pp.init.md),
which provides defaults of `seg.thresh = 0.05` and `miss.thresh = 0.1`.
These may be overridden through the `pars` argument.

## Setting markers aside

In the following, markers are pulled that are co-located with another
marker, that exhibit segregation distortion at a p-value below 0.02, or
that are missing more than 3% of their scores.

``` r

mapDHs <- pullCross(mapDH, type = "co.located")
mapDHs <- pullCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.02))
mapDHs <- pullCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.03))

names(mapDHs)
```

    [1] "geno"           "pheno"          "co.located"     "seg.distortion"
    [5] "missing"       

The object has acquired one element for each marker type. Each contains
a `table` summarising the markers concerned and a `data` matrix of their
scores, in the same genotype by marker layout as the linkage groups
themselves.

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

For the types `"seg.distortion"` and `"missing"`, the table combines
positional information with the output of
[`geno.table()`](https://rdrr.io/pkg/qtl/man/geno.table.html) from the
R/qtl package ([Broman and Wu 2014](#ref-br14)), so that the basis for
each marker’s removal may be inspected.

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

The table for `"co.located"` markers is of a different form. Markers are
assigned to bins, and within each bin the **first marker is retained in
the map as a reference**, the remainder being pulled. It is this
reference that permits the co-locating markers to be restored to their
correct position.

## Reconstruction

Reconstruction with `bychr = FALSE` combines all remaining markers and
rebuilds the map in its entirety, in the course of which the linkage
groups are renamed. This provides a reasonably stringent test of whether
the markers set aside can subsequently be restored to their correct
positions.

``` r

mapDHs <- mstmap(mapDHs, bychr = FALSE, dist.fun = "kosambi", trace = FALSE,
                 anchor = TRUE)
```

    Number of linkage groups: 24
    The size of the linkage groups are: 38  4   31  9   37  34  7   13  35  23  5   8   21  32  6   52  55  5   27  41  13  32  34  3   
    The number of bins in each linkage group: 38    4   31  9   37  34  7   13  35  23  5   8   21  32  6   52  55  5   27  41  13  32  34  3   

``` r

nmar(mapDHs)
```

     L.1 L.10 L.11 L.12 L.13 L.14 L.15 L.16 L.17 L.18 L.19  L.2 L.20 L.21 L.22 L.23 
      38   23    5    8   21   32    6   52   55    5   27    4   41   13   32   34 
    L.24  L.3  L.4  L.5  L.6  L.7  L.8  L.9 
       3   31    9   37   34    7   13   35 

## Reversal of the thresholds

The thresholds operate in opposite senses in the two functions, and this
is a frequent source of confusion.

> [`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
> removes a marker when its segregation distortion p-value lies
> **below** `seg.thresh`, whereas
> [`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
> restores it when that p-value lies **above** `seg.thresh`. The same
> applies to `miss.thresh`.

Pulling and pushing at an identical threshold therefore returns no
markers at all, and the function reports that there are none of that
type to push back. The two thresholds are intended to differ, the value
used when pushing being the less stringent of the two, so that a marker
is restored only if it is adequately behaved in the context of the
established map.

| Parameter | Pulled when | Restored when |
|----|----|----|
| `seg.thresh` | p-value **\<** threshold | p-value **\>** threshold |
| `miss.thresh` | missing proportion **\>** threshold | missing proportion **\<** threshold |

Markers were pulled above at `seg.thresh = 0.02` and
`miss.thresh = 0.03`, and are accordingly restored at 0.001 and 0.05
respectively.

``` r

mapDHs <- pushCross(mapDHs, type = "co.located")
mapDHs <- pushCross(mapDHs, type = "seg.distortion", pars = list(seg.thresh = 0.001))
mapDHs <- pushCross(mapDHs, type = "missing", pars = list(miss.thresh = 0.05))

names(mapDHs)
```

    [1] "geno"  "pheno"

Under these settings all markers of each type have been restored, and
the holding elements have consequently been removed from the object.

## Placement of the restored markers

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

Co-located markers are placed adjacent to their reference marker, as for
example `1D.m.4` beside `1D.m.3` on 1D. Markers held in the
`"seg.distortion"` and `"missing"` elements cannot be positioned so
precisely and are appended to the end of the linkage group to which they
are assigned; the distorted marker `6D.m.12` appears at the end of 6D.

> [`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
> assigns markers to linkage groups but does not order them. The map is
> not complete until MSTmap has been run again.

A final pass within each linkage group places all markers in their
optimal position. Setting `p.value = 2` suppresses clustering, so that
the composition of the groups cannot change.

``` r

mapDHs <- mstmap(mapDHs, bychr = TRUE, dist.fun = "kosambi", trace = FALSE,
                 anchor = TRUE, p.value = 2)
```

    Number of linkage groups: 1
    The size of the linkage groups are: 41  
    The number of bins in each linkage group: 41    
    Number of linkage groups: 1
    The size of the linkage groups are: 30  
    The number of bins in each linkage group: 29    
    Number of linkage groups: 1
    The size of the linkage groups are: 6   
    The number of bins in each linkage group: 5 
    Number of linkage groups: 1
    The size of the linkage groups are: 9   
    The number of bins in each linkage group: 8 
    Number of linkage groups: 1
    The size of the linkage groups are: 21  
    The number of bins in each linkage group: 21    
    Number of linkage groups: 1
    The size of the linkage groups are: 32  
    The number of bins in each linkage group: 32    
    Number of linkage groups: 1
    The size of the linkage groups are: 6   
    The number of bins in each linkage group: 6 
    Number of linkage groups: 1
    The size of the linkage groups are: 54  
    The number of bins in each linkage group: 53    
    Number of linkage groups: 1
    The size of the linkage groups are: 56  
    The number of bins in each linkage group: 55    
    Number of linkage groups: 1
    The size of the linkage groups are: 6   
    The number of bins in each linkage group: 5 
    Number of linkage groups: 1
    The size of the linkage groups are: 27  
    The number of bins in each linkage group: 27    
    Number of linkage groups: 1
    The size of the linkage groups are: 5   
    The number of bins in each linkage group: 4 
    Number of linkage groups: 1
    The size of the linkage groups are: 41  
    The number of bins in each linkage group: 41    
    Number of linkage groups: 1
    The size of the linkage groups are: 15  
    The number of bins in each linkage group: 14    
    Number of linkage groups: 1
    The size of the linkage groups are: 33  
    The number of bins in each linkage group: 32    
    Number of linkage groups: 1
    The size of the linkage groups are: 37  
    The number of bins in each linkage group: 36    
    Number of linkage groups: 1
    The size of the linkage groups are: 4   
    The number of bins in each linkage group: 3 
    Number of linkage groups: 1
    The size of the linkage groups are: 33  
    The number of bins in each linkage group: 33    
    Number of linkage groups: 1
    The size of the linkage groups are: 10  
    The number of bins in each linkage group: 9 
    Number of linkage groups: 1
    The size of the linkage groups are: 40  
    The number of bins in each linkage group: 40    
    Number of linkage groups: 1
    The size of the linkage groups are: 35  
    The number of bins in each linkage group: 34    
    Number of linkage groups: 1
    The size of the linkage groups are: 8   
    The number of bins in each linkage group: 7 
    Number of linkage groups: 1
    The size of the linkage groups are: 13  
    The number of bins in each linkage group: 13    
    Number of linkage groups: 1
    The size of the linkage groups are: 37  
    The number of bins in each linkage group: 36    

## Markers not previously present in the map

[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
admits a fourth type, `"unlinked"`, which distributes markers residing
in an unlinked group of the `geno` element into the established linkage
groups. This provides the mechanism by which new markers are introduced
into a completed map, and is described in [Worked example
II](https://drj001.github.io/ASMap/articles/worked-example-refinement.md)
in the context of fine mapping and of unknown linkage groups.

## Further reading

- [Diagnosing genotypes and
  markers](https://drj001.github.io/ASMap/articles/diagnostics.md) —
  determining which markers warrant setting aside
- [Worked example
  I](https://drj001.github.io/ASMap/articles/worked-example-construction.md)
  — the pull operation within a complete construction
- [`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md),
  [`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md),
  [`pp.init()`](https://drj001.github.io/ASMap/reference/pp.init.md) —
  full argument documentation

## References

Broman, K. W, and H Wu. 2014. *qtl: Tools for Analyzing QTL
Experiments*. <https://CRAN.R-project.org/package=qtl>.

Taylor, Julian, and David Butler. 2017. “R Package ASMap: Efficient
Genetic Linkage Map Construction and Diagnosis.” *Journal of Statistical
Software* 79 (6): 1–29. <https://doi.org/10.18637/jss.v079.i06>.
