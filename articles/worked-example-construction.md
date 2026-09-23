# Worked example I: construction

This chapter involves the complete linkage map construction process for
a barley Backcross population that contains 3024 markers genotyped on
326 individuals in an unconstructed marker set formatted as an R/qtl
object with class `"bc"`. The data is available from the R/ASMap package
by typing

``` r

data(mapBCu, package = "ASMap")
```

## Pre-construction

The construction of a linkage map does not usually just involve applying
a construction algorithm to a supplied set of genetic marker data. It is
always prudent to go through a pre-construction checklist to ensure that
the best quality genotypes/markers are being used to construct the
linkage map.

A non-exhaustive ordered checklist for an unconstructed marker set could
be

1.  Check missing allele scores across markers for each genotype as well
    as across genotypes for each marker. Markers or genotypes with a
    high proportion of missing information could indicate problems with
    the physical genotyping.

2.  Check for genetic clones or individuals that have a high proportion
    of matching allelic information between them.

3.  Check markers for excessive segregation distortion. Highly distorted
    markers may not map to unique locations.

4.  Check markers for switched alleles. These markers will not cluster
    or link well with other markers during the construction process and
    it is therefore preferred to repair their alignment before
    proceeding.

5.  Check for co-locating markers. For large linkage maps it would be
    more computationally efficient from a construction standpoint to
    temporarily omit markers that are co-located with other markers.

R/qtl provides a very simple graphical tool for checking the structure
of missing allele score across the genotypes and the markers.

``` r

plotMissing(mapBCu)
```

![Plot of the missing allele scores for the unconstructed map
mapBCu](worked-example-construction_files/figure-html/ex3-1.png)

Plot of the missing allele scores for the unconstructed map mapBCu

the figure below shows the resulting plot for this command. The darkest
lines on the plot indicate there are some genotypes with large amounts
of missing data. This could indicate poor physical genotyping of these
lines and they should be removed before proceeding. The plot also
reveals the markers have a large number of typed allele values across
the range of genotypes. The R/ASMap function
[`statGen()`](https://drj001.github.io/ASMap/reference/statGen.md) can
be used to identify the genotypes with a certain number of missing
values. These genotypes are then omitted using the usual functions
available in R/qtl.

``` r

sg <- statGen(mapBCu, bychr = FALSE, stat.type = "miss")
mapBC1 <- subset(mapBCu, ind = sg$miss < 1600)
```

From a map construction point of view, highly related individuals can
enhance segregation distortion of markers. It is therefore wise to
determine a course of action such as removal of individuals or the
creation of consensus genotypes before proceeding with any further
pre-construction diagnostics. The R/ASMap function
[`genClones()`](https://drj001.github.io/ASMap/reference/genClones.md)
discussed in section Genetic clones can be used to identify and report
genetic clones.

``` r

gc <- genClones(mapBC1, tol = 0.95)
gc$cgd
```

          G1    G2   coef match diff na.both na.one group
    1  BC045 BC039 0.9919  1466   12      97   1448     1
    2  BC052 BC039 1.0000  2423    0      98    502     1
    3  BC168 BC039 1.0000  2572    0      47    404     1
    4  BC052 BC045 0.9899  1476   15      94   1438     1
    5  BC168 BC045 0.9872  1620   21      44   1338     1
    6  BC168 BC052 1.0000  2577    0      36    410     1
    7  BC067 BC060 1.0000  2759    0       8    256     2
    8  BC135 BC086 1.0000  2743    0      17    263     3
    9  BC099 BC093 1.0000  2737    0      19    267     4
    10 BC120 BC117 1.0000  2699    0      22    302     5
    11 BC129 BC126 1.0000  2678    0      35    310     6
    12 BC204 BC138 1.0000  2691    0       6    326     7
    13 BC147 BC141 0.9996  2753    1      22    247     8
    14 BC193 BC144 1.0000  2771    0       7    245     9
    15 BC162 BC161 1.0000  2686    0      20    317    10
    16 BC205 BC190 1.0000  2886    0       7    130    11
    17 BC286 BC285 1.0000  2920    0       1    102    12
    18 BC325 BC314 1.0000  2911    0       4    108    13

The table shows 13 groups of genotypes that share a proportion of their
alleles greater than 0.95. The supplied additional statistics show that
the first group contains three pairs of genotypes that had matched pairs
of alleles from 1620 markers or less. These pairs also   1400 markers
where an allele was present for one genotype and missing for another.
Based on this, there is not enough evidence to suspect these pairs may
be clones and they are removed from the table. The
[`fixClones()`](https://drj001.github.io/ASMap/reference/fixClones.md)
function can then be used to form consensus genotypes for the remaining
groups of clones in the table.

``` r

cgd <- gc$cgd[-c(1,4,5),]
mapBC2 <- fixClones(mapBC1, cgd, consensus = TRUE)
levels(mapBC2$pheno[[1]])[grep("_", levels(mapBC2$pheno[[1]]))]
```

     [1] "BC039_BC052_BC168" "BC060_BC067"       "BC086_BC135"      
     [4] "BC093_BC099"       "BC117_BC120"       "BC126_BC129"      
     [7] "BC138_BC204"       "BC141_BC147"       "BC144_BC193"      
    [10] "BC161_BC162"       "BC190_BC205"       "BC285_BC286"      
    [13] "BC314_BC325"      

At this juncture it is wise to check the segregation distortion
statistics of the markers. Segregation distortion is phenomenon where
the observed allelic frequencies at a specific locus deviate from
expected allelic frequencies due to Mendelian genetics. It is well known
that this distortion can occur from physical laboratorial processes or
it may also occur in local genomic regions from underlying biological
and genetic mechanisms ([Lyttle 1991](#ref-lytt91)).

The level of segregation distortion, the allelic proportions and the
missing value proportion across the genome can be graphically
represented using the marker profiling function
[`profileMark()`](https://drj001.github.io/ASMap/reference/profileMark.md)
and the result is displayed in the figure below.

Setting `crit.val = "bonf"` annotates the markers that have a p-value
for the test of segregation distortion lower than the family wide
bonferroni adjusted alpha level of 0.05/no.of.markers on each of the
figures in each panel. Additional plotting parameters `layout = c(1,4)`,
`type="p"` and `cex = 0.5` are passed to the high level lattice plotting
function [`xyplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html) to
provide a more aesthetically pleasing plot. The plot indicates there are
numerous markers that are considered to be significantly distorted with
three highly distorted markers. The plot also indicates that the missing
value proportion of the markers does not exceed 20%.

``` r

profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"), crit.val = "bonf", layout = c(1,4), type = "l", cex = 0.5)
```

![For individual markers, the negative log10 p-value for the test of
segregation distortion, the proportion of each contributing allele and
the proportion of missing
values.](worked-example-construction_files/figure-html/ex8-1.png)

For individual markers, the negative log10 p-value for the test of
segregation distortion, the proportion of each contributing allele and
the proportion of missing values.

The highly distorted markers can easily be dropped using

``` r

mm <- statMark(mapBC2, stat.type = "marker")$marker$AB
mapBC3 <- drop.markers(mapBC2, c(markernames(mapBC2)[mm > 0.98],markernames(mapBC2)[mm < 0.2]))
```

Without a constructed map it is impossible to determine the origin of
the segregation distortion. However, the blind use of distorted markers
may also create linkage map construction problems. It may be more
sensible to place the distorted markers aside and construct the map with
less problematic markers. Once the linkage map is constructed the more
problematic markers can be introduced to determine whether they have a
useful or deleterious effect on the map. The R/ASMap functions
[`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
and
[`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
are designed to take advantage of this scenario. To showcase their use
in this example, they are also used to pull markers with 10-20% missing
values as well as co-located markers.

``` r

mapBC3 <- pullCross(mapBC3, type = "missing", pars = list(miss.thresh = 0.1))
mapBC3 <- pullCross(mapBC3, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
mapBC3 <- pullCross(mapBC3, type = "co.located")
names(mapBC3)
```

    [1] "geno"           "pheno"          "missing"        "seg.distortion"
    [5] "co.located"    

``` r

sum(ncol(mapBC3$missing$data),ncol(mapBC3$seg.dist$data),ncol(mapBC3$co.located$data))
```

    [1] 847

A total of 847 markers are removed and placed aside in their respective
elements and the map is now constructed with the remaining 2173 markers.

## MSTmap construction

The curated genetic marker data in `mapBC3` can now be constructed using
the `mstmap.cross` function available in R/ASMap.

``` r

mapBC4 <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12)
chrlen(mapBC4)
```

           L.1        L.2        L.3        L.4        L.5        L.6        L.7 
    304.910957 266.240647  78.982131 252.281760  33.226962 233.485952 153.315888 
           L.8        L.9 
    106.290403   6.657526 

By setting `bychr = FALSE` the complete set of marker data from `mapBC3`
is bulked and constructed from scratch. This construction involves the
clustering of markers to linkage groups and then the optimal ordering of
markers within each linkage group. An initial check of the figure below,
indicates for a population size of 309 the `p.value` should be set to
`1e-12` to ensure a 30cM threshold when clustering markers to linkage
groups. The newly constructed linkage map contains 9 linkage groups each
containing markers that are optimally ordered. The performance of the
MSTmap construction can be checked by plotting the heat map of pairwise
recombination fractions between markers and their pairwise LOD score of
linkage using the R/ASMap function
[`heatMap()`](https://drj001.github.io/ASMap/reference/heatMap.md)

``` r

heatMap(mapBC4, lmax = 70)
```

    Warning in heatMap(mapBC4, lmax = 70): Running est.rf.

![Heat map of the constructed linkage map
mapBC4.](worked-example-construction_files/figure-html/ex12-1.png)

Heat map of the constructed linkage map mapBC4.

As discussed in section Improved heat map, an aesthetic heat map is
attained when the heat on the upper triangle of the plot used for the
pairwise estimated recombination fractions matches the heat of the
pairwise LOD scores. the figure below displays the heat map and shows
this was achieved by setting `lmax = 70`. The heat map also shows
consistent heat across the markers within linkage groups indicating
strong linkage between nearby markers. The linkage groups appear to be
very distinctly clustered.

Although the heat map is indicating the construction process was
successful it does not highlight subtle problems that may be existing in
the constructed linkage map. As stated in section Genotype and
marker/interval profiling one of the key quality characteristics of a
well constructed linkage map is an appropriate recombination rate of the
the genotypes. For this barley Backcross population each line is
considered to be independent with an expected recombination rate of   14
across the genome. The current linkage map contains linkage groups that
exceed the theoretical cutoff of 200cM indicating there may be genotypes
with inflated recombination rates. This can easily be ascertained from
the
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md)
function available in R/ASMap.

``` r

pg <- profileGen(mapBC4, bychr = FALSE, stat.type = c("xo","dxo","miss"), id = "Genotype", xo.lambda = 14, layout = c(1,3), lty = 2, cex = 0.7)
```

![For individual genotypes, the number of recombinations, double
recombinations and missing values for
mapBC4](worked-example-construction_files/figure-html/ex15-1.png)

For individual genotypes, the number of recombinations, double
recombinations and missing values for mapBC4

the figure below show the number of recombinations, double recombination
and missing values for each of 309 genotypes. The plot also annotates
the genotypes that have recombination rates significantly above the
expected recombination rate of 14. A total of seven lines have
recombination rates above 20 and the plots also show that these lines
have excessive missing values. To ensure the extra list elements
`"co.located"`, `"seg.distortion"` and `"missing"` of the object are
subsetted and updated appropriately, the offending genotypes are removed
using the R/ASMap function
[`subsetCross()`](https://drj001.github.io/ASMap/reference/subsetCross.md).
The linkage map is then be reconstructed.

``` r

mapBC5 <- subsetCross(mapBC4, ind = !pg$xo.lambda)
mapBC6 <- mstmap(mapBC5, bychr = TRUE, dist.fun = "kosambi", trace = TRUE, p.value = 1e-12)
chrlen(mapBC6)
```

           L.1        L.2        L.3        L.4        L.5        L.6        L.7 
    229.578103 223.170605  65.806120 214.380402  27.297642 191.060846 140.114968 
           L.8        L.9 
     88.687053   4.875885 

Users can check the recombination rates of the remaining genotypes in
the re-constructed map are now within respectable limits. As a result
the lengths of the linkage groups have dropped dramatically.

It is also useful to graphically display statistics of the markers and
intervals of the current constructed linkage map. For example, the
figure below shows the marker profiles of the -log10 p-value for the
test of segregation distortion, the allele proportions and the number of
double crossovers. It also displays the interval profile of the number
of recombinations occurring between adjacent markers. This plot reveals
many things that are useful for the next phase of the construction
process.

``` r

profileMark(mapBC6, stat.type = c("seg.dist","prop","dxo","recomb"), layout = c(1,5), type = "l")
```

![Marker profiles of the -log10 p-value for the test of segregation
distortion, allele proportions and the number of double of crossovers as
well as the interval profile of the number of recombinations between
adjacent markers in
mapBC6.](worked-example-construction_files/figure-html/ex18-1.png)

Marker profiles of the -log10 p-value for the test of segregation
distortion, allele proportions and the number of double of crossovers as
well as the interval profile of the number of recombinations between
adjacent markers in mapBC6.

The plot instantly reveals the success of the map construction process
with no more than one double crossover being found at any marker and
very few being found in total. The plot also reveals the extent of the
biological distortion that can occur within a linkage group. A close
look at the segregation distortion and allele proportion plots shows the
linkage group L.3 and the short linkage group L.5 have profiles that
could be joined if the linkage groups were merged. In addition, L.8 and
L.9 also have profiles that could be joined if the linkage groups were
combined. This will be discussed in more detail in the next section.

Lyttle, T. W. 1991. “Segregation Distorters.” *Annula Reviews of
Genetics* 25: 511–86.
