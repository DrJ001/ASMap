# Constructing a linkage map

This chapter explores the R/ASMap functions in greater depth and shows
how they can be used to efficiently explore, manipulate and construct
genetic linkage maps. It will also showcase the graphical tools that
will allow efficient diagnosis of linkage map problems post
construction.

The package contains multiple data sets listed as follows

- **mapDH**:

  : A constructed genetic linkage map for a Doubled Haploid population
  in the form of an R/qtl object. The genetic linkage map contains a
  total of 599 markers spanning 23 linkage groups genotyped across 218
  individuals. The linkage map contains a small set of co-located
  markers and a small set of markers with excessive segregation
  distortion

- **mapDHf**:

  : An unconstructed version of `mapDH` in the form a data frame. The
  data frame has dimensions 599 $`\times`$ 218 and the rows (markers)
  have been randomized.

- **mapBCu**:

  : An unconstructed set of markers for a backcross population in the
  form of an R/qtl object. The marker set contains a total of 3023
  markers genotyped on 326 individuals. This marker set can be used in
  conjunction with the detailed process presented in chapter to
  construct a genetic linkage map.

- **mapBC**:

  : A constructed linkage map of `mapBCu` in the form of an R/qtl
  object. The linkage map contains a total of 3019 markers genotyped on
  300 individuals.

- **mapF2**:

  : A simulated linkage map for a self-pollinated F2 population
  consisting of 700 markers spanning 7 linkage groups genotyped across
  250 individuals.

Each of these data sets is accessible using commands similar to

``` r

data(mapDHf, package = "ASMap")
data(mapDH, package = "ASMap")
data(mapBCu, package = "ASMap")
```

## Map construction functions

The R/ASMap package contains two linkage map construction functions that
allow users to fully utilize the MSTmap parameters listed at
[http://alumni.cs.ucr.edu/~yonghui/mstmap.html](http://alumni.cs.ucr.edu/~yonghui/mstmap.md).
Some additional information on aspects of the MSTmap algorithm and the
appropriate use of the arguments `mvest.bc` and `detectBadData` is given
in Chapter .

### `mstmap.data.frame()`

The first of these functions allows users to input a data frame of
genetic markers ready for construction. For a more detailed explanation
of the arguments users should consult the help documentation found by
typing
[`?mstmap.data.frame`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
in R.

``` r

mstmap.data.frame(object, pop.type = "DH", dist.fun = "kosambi",
      objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
      noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData = FALSE,
      as.cross = TRUE, return.imputed = TRUE, trace = FALSE, ...)
```

The explicit form of the data frame `object` is borne from the syntax of
the marker file required for using the stand alone MSTmap software. It
must have markers in rows and genotypes in columns. Marker names are
required to be in the `rownames` component of the object with genotype
names residing in the `names`. Spaces in any of the marker or genotype
names should be avoided but will be replaced with a “-” if found. Each
of the columns of the data frame must be of class `"character"` (not
factors). If converting from a matrix, this can easily be achieved by
using the `stringAsFactors = FALSE` argument for any `data.frame`
method.

The available populations that can be passed to `pop.type` are `"BC"`
Backcross, `"DH"` Doubled Haploid, `"ARIL"` Advanced Recombinant Inbred
and `"RILn"` Recombinant Inbred with n levels of selfing. The allelic
content of the markers in the `object` must be explicitly adhered to.
For `pop.type` `"BC"`, `"DH"` or `"ARIL"` the two allele types should be
represented as (`"A"` or `"a"`) and (`"B"` or `"b"`). Thus for
`pop.type = "ARIL"` it is assumed the minimal number of heterozygotes
have been set to missing. For non-advanced RIL populations
(`pop.type = "RILn"`) phase unknown heterozygotes should be represented
as `"X"`. For all populations, missing marker scores should be
represented as (`"U"` or `"-"`).

Users need to be aware that the `p.value` argument plays an important
role in determining the clustering of markers to distinct linkage
groups. Section Clustering shows the separation of marker groups is
highly dependent on the the number of individuals in the population. For
this reason, some trial and error may be required to determine an
appropriate `p.value` for the linkage map being constructed.

Although this function contains arguments that utilize the complete set
of available MSTmap parameters it is less flexible than its sister
function
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
(see section textttmstmap.cross()) that uses the flexible structure of
an R/qtl `"cross"` object. For this reason, it is recommended that users
set `as.cross = TRUE` to ensure the constructed object is returned as a
R/qtl cross object with an appropriate class structure. For population
types `"BC"` and `"DH"` the class of the constructed object is given
`"bc"` and `"dh"` respectively. For `"RILn"` the **qtl** package
conversion function `convert2bcsft` is used to ensure the class of the
object is assigned `"bcsft"` with arguments `F.gen = n` and
`BC.gen = 0`. For `"ARIL"` populations the constructed object is given
the class `"riself"`. The correct assignation of these classes ensures
the objects can be used synergistically with the suite of functions
available in the R/qtl package as well as other functions in the R/ASMap
package.

The R/ASMap package contains an unconstructed Doubled Haploid marker set
`mapDHf` with 599 markers genotyped across 218 individuals. The marker
set is formatted correctly for input into the
[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
function.

``` r

testd <- mstmap(mapDHf, dist.fun = "kosambi", trace = TRUE, as.cross = TRUE)
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

As `as.cross = TRUE` the usual functions available in the R/qtl package
are available for use on the returned object.

### `mstmap.cross()`

The second linkage map construction function allows users to input an
unconstructed or constructed linkage map in the form of an R/qtl cross
object. See
[`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
for a more detailed description.

``` r

mstmap.cross(object, chr, id = "Genotype", bychr = TRUE,
       suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
       objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
       noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE, detectBadData =
       FALSE, return.imputed = FALSE, trace = FALSE, ...)
```

The cross `object` needs to inherit from one of the allowable classes
available in the R/qtl package, namely `"bc","dh","riself","bcsft"`
where `"bc"` is a Backcross `"dh"` is Doubled Haploid, `"riself"` is an
advanced Recombinant Inbred and `"bcsft"` is a Backcross/Self.

It is important to understand how these classes are encoded into the
object for the specific populations. If
[`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html) is used to
read in any bi-parental populations it will be given the class `"bc"`.
Doubled Haploid populations can be changed to `"dh"` just by changing
the class. For the purpose of linkage map construction, both classes
`"bc"` and `"dh"` will produce equivalent results. For non-advanced
Recombinant Inbred populations markers are required to be fully
informative (i.e. contain 3 distinct allele types such as AA, BB for
parental homozygotes and AB for phase unknown heterozygotes) and the use
of [`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html) will
result in the cross object being given a class `"f2"`. The level of
selfing is required to be encoded into the object by applying one of the
two conversion functions available in the R/qtl package. For a
population that has been generated by selfing $`n`$ times, the
conversion function `convertbcsft` can be used by setting the arguments
`F.gen = n` and `BC.gen = 0`. This will attach a class `"bcsft"` to the
object. Populations that are genuine advanced RILs can be converted
using the `convert2riself` function. This function will replace any
remaining heterozygosity, if it exists, with missing values and attach
the class `"riself"` to the object.

Similar to the
[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
function, users need to be aware that the `p.value` argument is highly
dependent on the number of individuals in the population and may require
some trial and error to ascertain an appropriate value. After
construction the cross object is returned with an identical class
structure as the inputted object. All R/qtl and R/ASMap package
functions can be used synergistically with this object.

#### Examples

The constructed linkage map `mapDH` available in the R/ASMap package
will be used to showcase the flexibility of this function. Before
attempting re-construction, some preliminary output of `mapDH` is
presented.

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

The output shows that there are 23 groups that have been appropriately
been assigned linkage group or chromosome names. The markers within each
linkage group have been named according to the order of the markers.

**Example 1: Completely construct or reconstruct a linkage map.**

To completely re-construct this map set the argument `bychr = FALSE`.
This will bulk the genetic data from all linkage groups, re-cluster the
markers into groups and then optimally order the markers within each
linkage group. This linkage map is small so these two tasks happen
almost instantaneously.

``` r

mapDHa <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi", trace = TRUE)
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

The reconstructed map contains 24 linkage groups with the extra linkage
group coming from a minor split in 4A. As the linkage map is constructed
from scratch, it assumed that former linkage group names are no longer
required. A standard “L.” prefix is provided for the new linkage group
names. This example also indicates that MSTmap by default does not
respect the inputted marker order.

**Example 2: Optimally order markers within linkage groups of an
established map.**

It may be necessary to only perform marker ordering within already
established linkage groups. This can be achieved by setting
`bychr = TRUE`. In some cases it also be preferable to ensure the marker
orders of the linkage groups are respected and this can be achieved by
setting `anchor = TRUE.`

``` r

mapDHb <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE, trace = TRUE)
nmar(mapDHb)
```

      1A  1B1  1B2   1D   2A   2B  2D1  2D2   3A   3B   3D 4A.1 4A.2   4B   4D   5A 
      41    5   33   10   40   35    8   13   37   30    6    9   21   32    6   54 
      5B   5D   6A   6B   6D   7A   7B   7D 
      56    6   27   41   15   33   37    4 

This map is identical to `mapDH` with the exception that chromosome 4A
has been split into two linkage groups. As `bychr = TRUE` the function
understands the origin of the linkage group was 4A and consequently uses
it as a prefix in the naming of the two new linkage groups.

**Example 3: Optimally order markers within linkage groups of an
established map without breaking linkage groups.**

The splitting of the linkage groups in the last example only occurred
due to choice of default `p.value = 1e-06` set in the function. A slight
change to this `p.value` will ensure that 4A remains linked during the
algorithm.

``` r

mapDHc <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", anchor = TRUE, trace = TRUE, p.value = 1e-04)
nmar(mapDHc)
```

     1A 1B1 1B2  1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D 
     41   5  33  10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15 
     7A  7B  7D 
     33  37   4 

An identical result can be achieved by setting `p.value = 2` or any
number greater than 1. Doing this instructs the MSTmap algorithm to not
split linkage groups regardless of how weak the linkages are between
markers within any group. Users should be aware that if this latter
method is used then linkage groups that contain groups of markers
separated by a substantial distance (i.e. very weak linkages) may suffer
from local orientation problems.

**Example 4: Reconstruct map within predefined linkage groups of an
established map.**

There may only be a need to reconstruct a predefined set of linkage
groups. By setting the `chr` argument and `bychr = FALSE` users can
determine which linkage groups require the complete reconstruction using
MSTmap.

``` r

mapDHd <- mstmap(mapDH, chr = names(mapDH$geno)[1:3], bychr = FALSE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-04)
nmar(mapDHd)
```

     1D  2A  2B 2D1 2D2  3A  3B  3D  4A  4B  4D  5A  5B  5D  6A  6B  6D  7A  7B  7D 
     10  40  35   8  13  37  30   6  30  32   6  54  56   6  27  41  15  33  37   4 
    L.1 L.2 L.3 
     41   5  33 

Again, the algorithm assumes that the original linkage group names are
no longer necessary and defines new ones. An obvious extension of this
example is to set `bychr = TRUE` and then the algorithm will order
markers within the predefined linkage groups stipulated by `chr`.
Similar to the previous example, an appropriate choice of `p.value` will
ensure that linkage groups remain unbroken.
