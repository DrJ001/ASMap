# Constructing a linkage map

ASMap provides two construction functions, both methods of the
[`mstmap()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
generic. They differ only in the object they accept: a data frame of
marker scores, or an R/qtl `cross` object. Both expose the full set of
MSTmap parameters documented at
[http://alumni.cs.ucr.edu/~yonghui/mstmap.html](http://alumni.cs.ucr.edu/~yonghui/mstmap.md).

The algorithm underlying them is described in [How the MSTmap algorithm
works](https://drj001.github.io/ASMap/articles/mstmap-algorithm.md), and
the practical behaviour of the `mvest.bc` and `detectBadData` arguments
in [Notes on the MSTmap
algorithm](https://drj001.github.io/ASMap/articles/algorithm-notes.md).

## Example data

The package supplies five data sets, described in full on their
reference pages. The datasets are not lazy-loaded and require an
explicit [`data()`](https://rdrr.io/r/utils/data.html) call.

| Data set | Form | Individuals | Markers | Description |
|----|----|----|----|----|
| [`mapDH`](https://drj001.github.io/ASMap/reference/mapDH.md) | cross | 218 | 599 | Constructed doubled haploid map, 23 linkage groups, containing a few co-located and distorted markers |
| [`mapDHf`](https://drj001.github.io/ASMap/reference/mapDHf.md) | data frame | 218 | 599 | Unconstructed form of `mapDH`, with the marker rows randomised |
| [`mapBCu`](https://drj001.github.io/ASMap/reference/mapBCu.md) | cross | 326 | 3023 | Unconstructed backcross marker set |
| [`mapBC`](https://drj001.github.io/ASMap/reference/mapBC.md) | cross | 300 | 3019 | Constructed form of `mapBCu` |
| [`mapF2`](https://drj001.github.io/ASMap/reference/mapF2.md) | cross | 250 | 700 | Simulated selfed F2 map, 7 linkage groups |

``` r

data(mapDHf, package = "ASMap")
data(mapDH, package = "ASMap")
data(mapBCu, package = "ASMap")
```

## Construction from a data frame

[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
accepts a data frame of marker scores. Its required layout derives from
the marker file used by the standalone MSTmap software: **markers in
rows and genotypes in columns**, with marker names in the `rownames` and
genotype names in the `names`. Every column must be of class
`"character"` rather than a factor, which is most easily ensured through
the `stringsAsFactors = FALSE` argument of any `data.frame` method.
Spaces in marker or genotype names should be avoided; where found they
are replaced by a hyphen.

The allelic coding is strict, and depends on the population type given
to `pop.type`.

| `pop.type` | Population | Allele coding |
|----|----|----|
| `"BC"` | Backcross | `"A"` or `"a"`, and `"B"` or `"b"` |
| `"DH"` | Doubled haploid | `"A"` or `"a"`, and `"B"` or `"b"` |
| `"ARIL"` | Advanced recombinant inbred | `"A"` or `"a"`, and `"B"` or `"b"`; heterozygotes assumed already set to missing |
| `"RILn"` | Recombinant inbred, *n* levels of selfing | As above, with phase-unknown heterozygotes coded `"X"` |

Missing scores are denoted `"U"` or `"-"` for all population types.

`mapDHf` is an unconstructed doubled haploid marker set of 599 markers
scored on 218 individuals, formatted for input to this function.

``` r

testd <- mstmap(mapDHf, dist.fun = "kosambi", trace = FALSE, as.cross = TRUE)
```

    Number of linkage groups: 24
    The size of the linkage groups are: 56  54  35  41  4   37  6   30  40  27  37  13  15  10  21  6   32  33  33  41  6   9   5   8   
    The number of bins in each linkage group: 55    53  34  41  3   36  6   29  40  27  36  13  14  9   21  5   32  33  32  41  5   8   4   7   

``` r

nmar(testd)
```

     L1  L2  L3  L4  L5  L6  L7  L8  L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 
     56  54  35  41   4  37   6  30  40  27  37  13  15  10  21   6  32  33  33  41 
    L21 L22 L23 L24 
      6   9   5   8 

``` r

chrlen(testd)
```

            L1         L2         L3         L4         L5         L6         L7 
    142.806369 181.935011  97.376797 115.405242  21.660902 108.162200  13.256931 
            L8         L9        L10        L11        L12        L13        L14 
    133.885978 153.235513  98.705941 145.839837  60.441336  72.211612  81.224002 
           L15        L16        L17        L18        L19        L20        L21 
     98.211685   8.738474 109.472813  60.913640 134.491476 103.954898  18.019891 
           L22        L23        L24 
      7.801089   7.810474  19.196583 

Setting `as.cross = TRUE` is recommended, since it returns an R/qtl
cross object of the appropriate class and thereby makes the facilities
of both packages available for subsequent work. The class assigned
depends on the population type.

| `pop.type` | Class assigned | Mechanism |
|----|----|----|
| `"BC"` | `"bc"` | Direct |
| `"DH"` | `"dh"` | Direct |
| `"ARIL"` | `"riself"` | Direct |
| `"RILn"` | `"bcsft"` | Via [`qtl::convert2bcsft()`](https://rdrr.io/pkg/qtl/man/qtl-internal.html) with `F.gen = n`, `BC.gen = 0` |

With `as.cross = FALSE` the map is returned instead as a data frame
carrying additional columns for the linkage group, marker position and
genetic distance.

## Construction from a cross object

[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
accepts an unconstructed or constructed map held in an R/qtl cross
object, and is the more flexible of the two. The object must inherit one
of the classes `"bc"`, `"dh"`, `"riself"` or `"bcsft"`, a restriction
that guards against attempting construction for the more complex
population types R/qtl supports.

How these classes arise matters in practice.
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html) assigns
the class `"bc"` to any bi-parental population. Doubled haploid
populations may be assigned `"dh"` simply by changing the class; for the
purpose of construction the two are equivalent.

Non-advanced recombinant inbred populations require fully informative
markers, that is three distinct allele types, and
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html) will
assign them the class `"f2"`. The level of selfing must then be encoded
by conversion:

| Population | Conversion |
|----|----|
| Selfed $`n`$ times | `qtl::convert2bcsft(F.gen = n, BC.gen = 0)`, giving class `"bcsft"` |
| Genuine advanced RIL | [`qtl::convert2riself()`](https://rdrr.io/pkg/qtl/man/convert2riself.html), giving class `"riself"` and replacing any residual heterozygosity with missing values |

The object is returned with the same class structure it was given, so
that all functions of both packages remain applicable.

> For both functions, the `p.value` argument governs the separation of
> markers into linkage groups and is strongly dependent on population
> size. Some experimentation is generally required. See [Clustering
> markers into linkage
> groups](https://drj001.github.io/ASMap/articles/mstmap-algorithm.md)
> for the relationship, which may be examined directly with
> [`pValue()`](https://drj001.github.io/ASMap/reference/pValue.md).

## Modes of reconstruction

The constructed map `mapDH` illustrates the flexibility of
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md).
Its 23 linkage groups are named, and the markers within each are named
according to their order.

``` r

nmar(mapDH)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D 
     41   5  33  10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15 
     7A  7B  7D 
     33  37   4 

``` r

pull.map(mapDH)[[4]]
```

       1D.m.1    1D.m.2    1D.m.3    1D.m.4    1D.m.5    1D.m.6    1D.m.7    1D.m.8 
     0.000000  2.285816  2.742698  2.742698 21.465574 25.652893 57.200980 58.104090 
       1D.m.9   1D.m.10 
    58.568087 81.061896 
    attr(,"class")
    [1] "A"

The four modes below are selected through `bychr`, `chr` and `anchor`.

| Requirement                     | `bychr` | `chr` | Effect on linkage groups      |
|---------------------------------|---------|-------|-------------------------------|
| Rebuild the entire map          | `FALSE` | unset | Re-clustered and renamed      |
| Re-order within existing groups | `TRUE`  | unset | Retained, possibly split      |
| Rebuild nominated groups only   | `FALSE` | set   | Nominated groups re-clustered |
| Re-order nominated groups only  | `TRUE`  | set   | Nominated groups retained     |

### Complete reconstruction

Setting `bychr = FALSE` combines the marker data from all linkage
groups, re-clusters the markers and orders them within each resulting
group.

``` r

mapDHa <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi", trace = FALSE)
```

    Number of linkage groups: 24
    The size of the linkage groups are: 41  5   33  10  40  35  8   13  37  30  6   9   21  32  6   54  56  6   27  41  15  33  37  4   
    The number of bins in each linkage group: 41    4   33  9   40  34  7   13  36  29  5   8   21  32  6   53  55  5   27  41  14  32  36  3   

``` r

nmar(mapDHa)
```

     L.1 L.10 L.11 L.12 L.13 L.14 L.15 L.16 L.17 L.18 L.19  L.2 L.20 L.21 L.22 L.23 
      41   30    6    9   21   32    6   54   56    6   27    5   41   15   33   37 
    L.24  L.3  L.4  L.5  L.6  L.7  L.8  L.9 
       4   33   10   40   35    8   13   37 

``` r

pull.map(mapDHa)[[4]]
```

      4A.m.9   4A.m.7   4A.m.8   4A.m.6   4A.m.5   4A.m.3   4A.m.4   4A.m.2 
    0.000000 1.835687 1.835687 2.294415 2.753144 3.670678 4.129406 6.424595 
      4A.m.1 
    7.801089 
    attr(,"class")
    [1] "A"

The reconstructed map contains 24 linkage groups, the additional group
arising from a minor split in 4A. Because the map is built from first
principles the former linkage group names are discarded and a standard
`"L."` prefix applied. This example also shows that MSTmap does not, by
default, respect the marker order it was given.

### Ordering within established groups

Where only the order within existing linkage groups is at issue, set
`bychr = TRUE`. Setting `anchor = TRUE` additionally causes the inputted
marker orders to be respected.

``` r

mapDHb <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE, trace = FALSE)
```

    Number of linkage groups: 1
    The size of the linkage groups are: 41  
    The number of bins in each linkage group: 41    
    Number of linkage groups: 1
    The size of the linkage groups are: 5   
    The number of bins in each linkage group: 4 
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
    Number of linkage groups: 1
    The size of the linkage groups are: 30  
    The number of bins in each linkage group: 29    
    Number of linkage groups: 1
    The size of the linkage groups are: 6   
    The number of bins in each linkage group: 5 
    Number of linkage groups: 2
    The size of the linkage groups are: 9   21  
    The number of bins in each linkage group: 8 21  
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

``` r

nmar(mapDHb)
```

      1A  1B1  1B2   1D   2A   2B  2D1  2D2   3A   3B   3D 4A.1 4A.2   4B   4D   5A 
      41    5   33   10   40   35    8   13   37   30    6    9   21   32    6   54 
      5B   5D   6A   6B   6D   7A   7B   7D 
      56    6   27   41   15   33   37    4 

The result is identical to `mapDH` except that 4A has been split in two.
As `bychr = TRUE`, the function recognises the origin of the group and
uses 4A as a prefix in naming the two new groups.

### Ordering without splitting

The split above arose solely from the default `p.value = 1e-06`.
Relaxing it keeps 4A intact.

``` r

mapDHc <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE,
                 trace = FALSE, p.value = 1e-04)
```

    Number of linkage groups: 1
    The size of the linkage groups are: 41  
    The number of bins in each linkage group: 41    
    Number of linkage groups: 1
    The size of the linkage groups are: 5   
    The number of bins in each linkage group: 4 
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
    Number of linkage groups: 1
    The size of the linkage groups are: 30  
    The number of bins in each linkage group: 29    
    Number of linkage groups: 1
    The size of the linkage groups are: 6   
    The number of bins in each linkage group: 5 
    Number of linkage groups: 1
    The size of the linkage groups are: 30  
    The number of bins in each linkage group: 29    
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

``` r

nmar(mapDHc)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D 
     41   5  33  10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15 
     7A  7B  7D 
     33  37   4 

The same result follows from `p.value = 2`, or any value greater than
unity, which instructs the algorithm not to split linkage groups however
weak the linkages within them.

> Suppressing splitting entirely carries a cost. A linkage group
> containing sets of markers separated by a substantial distance, and
> therefore weakly linked, may suffer local orientation problems.

### Reconstruction of nominated groups

Supplying `chr` together with `bychr = FALSE` restricts complete
reconstruction to the nominated linkage groups.

``` r

mapDHd <- mstmap(mapDH, chr = names(mapDH$geno)[1:3], bychr = FALSE,
                 dist.fun = "kosambi", trace = FALSE, p.value = 1e-04)
```

    Number of linkage groups: 3
    The size of the linkage groups are: 41  5   33  
    The number of bins in each linkage group: 41    4   33  

``` r

nmar(mapDHd)
```

     1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D  7A  7B  7D 
     10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15  33  37   4 
    L.1 L.2 L.3 
     41   5  33 

Here too the original group names are discarded and new ones defined.
Setting `bychr = TRUE` instead orders markers within the groups
nominated by `chr`, and an appropriate `p.value` will keep those groups
unbroken.

## Further reading

- [How the MSTmap algorithm
  works](https://drj001.github.io/ASMap/articles/mstmap-algorithm.md) —
  the basis of `objective.fun`, `p.value`, `mvest.bc` and
  `detectBadData`
- [Pulling and pushing
  markers](https://drj001.github.io/ASMap/articles/pulling-and-pushing.md)
  — curating the marker set before construction
- [Worked example
  I](https://drj001.github.io/ASMap/articles/worked-example-construction.md)
  — construction applied to an unconstructed barley backcross
- [`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md),
  [`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
  — full argument documentation

## References
