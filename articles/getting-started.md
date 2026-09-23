# Getting started with ASMap

This article takes a set of raw marker scores and builds a finished
linkage map from them, in about five minutes. Along the way it covers
the one argument you will have to think about, and how to check that the
result is any good.

ASMap does this with the MSTmap algorithm, which clusters markers into
linkage groups and orders them within each group in a single pass — no
separate “rippling” step, and no overnight wait.

## Installation

``` r

install.packages("ASMap")

# development version
# install.packages("devtools")
devtools::install_github("DrJ001/ASMap")
```

## Building a map from marker scores

`mapDHf` is an unconstructed doubled haploid wheat marker set, in the
layout the construction function expects: **markers in rows, genotypes
in columns**, with marker names in the `rownames` and genotype names in
the `names`.

> The package datasets are not lazy-loaded, so every example needs an
> explicit [`data()`](https://rdrr.io/r/utils/data.html) call.
> Forgetting it gives `object 'mapDHf' not found`.

``` r

library(ASMap)
data(mapDHf)

dim(mapDHf)
```

    [1] 599 218

``` r

mapDHf[1:5, 1:6]
```

            DH1 DH2 DH3 DH4 DH5 DH6
    5B.m.2    B   B   B   A   A   B
    5A.m.52   B   B   A   A   B   A
    2B.m.1    A   B   A   B   B   A
    6B.m.3    B   B   A   B   A   A
    5A.m.43   B   B   A   B   B   A

One call clusters the markers into linkage groups, orders them within
each group and estimates the genetic distances.

``` r

map <- mstmap(mapDHf, pop.type = "DH", dist.fun = "kosambi", p.value = 1e-06)
```

    Number of linkage groups: 24
    The size of the linkage groups are: 56  54  35  41  4   37  6   30  40  27  37  13  15  10  21  6   32  33  33  41  6   9   5   8   
    The number of bins in each linkage group: 55    53  34  41  3   36  6   29  40  27  36  13  14  9   21  5   32  33  32  41  5   8   4   7   

``` r

nchr(map)
```

    [1] 24

``` r

round(chrlen(map), 1)
```

       L1    L2    L3    L4    L5    L6    L7    L8    L9   L10   L11   L12   L13 
    142.8 181.9  97.4 115.4  21.7 108.2  13.3 133.9 153.2  98.7 145.8  60.4  72.2 
      L14   L15   L16   L17   L18   L19   L20   L21   L22   L23   L24 
     81.2  98.2   8.7 109.5  60.9 134.5 104.0  18.0   7.8   7.8  19.2 

The result is an R/qtl `cross` object, so the whole of
[qtl](https://CRAN.R-project.org/package=qtl) is available to you
immediately.

``` r

summary(map)
```

    Warning in summary.cross(map): Some markers at the same position on chr
    L1,L2,L3,L5,L6,L8,L11,L13,L14,L16,L19,L21,L22,L23,L24; use jittermap().

        Doubled haploids

        No. individuals:    218 

        No. phenotypes:     1 
        Percent phenotyped: 100 

        No. chromosomes:    24 
            Autosomes:      L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 
                            L17 L18 L19 L20 L21 L22 L23 L24 

        Total markers:      599 
        No. markers:        56 54 35 41 4 37 6 30 40 27 37 13 15 10 21 6 32 33 33 41 
                            6 9 5 8 
        Percent genotyped:  99.8 
        Genotypes (%):      AA:51.2  BB:48.8 

``` r

plotMap(map)
```

![](getting-started_files/figure-html/qtl-1.png)

## Starting from an existing cross object

If the markers are already in a `cross` object, the same generic accepts
it. `bychr = FALSE` bulks every marker and rebuilds the map from
scratch, ignoring the existing linkage groups.

``` r

data(mapDH)
map2 <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi", p.value = 1e-06)
```

    Number of linkage groups: 24
    The size of the linkage groups are: 41  5   33  10  40  35  8   13  37  30  6   9   21  32  6   54  56  6   27  41  15  33  37  4   
    The number of bins in each linkage group: 41    4   33  9   40  34  7   13  36  29  5   8   21  32  6   53  55  5   27  41  14  32  36  3   

``` r

nchr(map2)
```

    [1] 24

Setting `bychr = TRUE` instead keeps the existing groups and only
re-orders the markers within them. Two combinations are worth committing
to memory:

``` r

# cluster and order from scratch, ignoring existing linkage groups
mstmap(cross, bychr = FALSE, p.value = 1e-12)

# re-order within existing linkage groups, never splitting them
mstmap(cross, bychr = TRUE, p.value = 2)
```

A `p.value` greater than 1 switches clustering off altogether, which is
what makes the second idiom work.

## Choosing `p.value`

`p.value` sets the threshold at which markers are split into separate
linkage groups, and the right value depends strongly on population size
— larger populations need a smaller one. Expect to try a few.
[`pValue()`](https://drj001.github.io/ASMap/reference/pValue.md) plots
the relationship, so the choice can be made deliberately rather than by
trial and error.

``` r

pValue(dist = c(25, 30, 35, 40), pop.size = 100:500, map.function = "kosambi")
```

![](getting-started_files/figure-html/pvalue-1.png)

For a population of around 300, a `p.value` of `1e-12` gives a threshold
of roughly 30 cM before markers are linked into the same group.

## Checking the result

[`heatMap()`](https://drj001.github.io/ASMap/reference/heatMap.md) shows
pairwise recombination fractions and LOD scores of linkage together,
each with its own legend. Consistent heat within a linkage group means
the markers near each other really are linked.

``` r

heatMap(mapDH, lmax = 50)
```

![](getting-started_files/figure-html/heat-1.png)

## Where to go next

| Article | What it covers |
|----|----|
| [How the MSTmap algorithm works](https://drj001.github.io/ASMap/articles/mstmap-algorithm.md) | The minimum spanning tree formulation, clustering and marker ordering |
| [Constructing a linkage map](https://drj001.github.io/ASMap/articles/constructing-a-map.md) | Both construction functions and every MSTmap parameter |
| [Pulling and pushing markers](https://drj001.github.io/ASMap/articles/pulling-and-pushing.md) | Setting problematic markers aside and reintroducing them |
| [Diagnosing genotypes and markers](https://drj001.github.io/ASMap/articles/diagnostics.md) | Profiling individuals, markers and intervals; genetic clones |
| [Heat maps](https://drj001.github.io/ASMap/articles/heat-maps.md) | Reading and tuning the recombination fraction and LOD display |
| [Manipulating linkage maps](https://drj001.github.io/ASMap/articles/manipulating-maps.md) | Breaking, merging, subsetting and combining maps |
| [Worked example I: construction](https://drj001.github.io/ASMap/articles/worked-example-construction.md) | A real barley backcross, from raw markers to a constructed map |
| [Worked example II: refinement](https://drj001.github.io/ASMap/articles/worked-example-refinement.md) | Pushing markers back, merging groups, fine mapping |
| [Notes on the MSTmap algorithm](https://drj001.github.io/ASMap/articles/algorithm-notes.md) | Distance calculations, `mvest.bc` and `detectBadData` |

If you use the package in published work, please cite it —
`citation("ASMap")`, or Taylor and Butler (2017), *Journal of
Statistical Software* **79**(6), 1–29,
[doi:10.18637/jss.v079.i06](https://doi.org/10.18637/jss.v079.i06).
