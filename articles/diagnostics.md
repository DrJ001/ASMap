# Diagnosing genotypes and markers

## Genotype and marker/interval profiling

To provide a complete system for efficient linkage map construction and
diagnosis the R/ASMap package contains functions that calculate linkage
map statistics across the markers for every genotype as well as across
the markers/intervals of the genome. These statistics can then be
profiled across simultaneous R/lattice panels for quick diagnosis of
linkage map attributes.

``` r

statGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"), id = "Genotype")
profileGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = NULL, ...)
statMark(cross, chr, stat.type = c("marker", "interval"), map.function = "kosambi")
profileMark(cross, chr, stat.type = "marker", use.dist = TRUE, map.function = "kosambi", crit.val = NULL, display.markers = FALSE, mark.line = FALSE, ...)
```

The [`statGen()`](https://drj001.github.io/ASMap/reference/statGen.md)
and
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md)
are functions for calculation and profiling of the statistics for the
genotypes across the chosen marker set defined by `chr`. Specifically,
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md)
actually calls
[`statGen()`](https://drj001.github.io/ASMap/reference/statGen.md) to
obtain the statistics for profiling and also returns the statistics
invisibly after plotting. Setting `bychr = TRUE` will also ensure that
the genotype profiles are plotted individually for each linkage groups
given `chr`. The current statistics that can be calculated for each
genotype include

$`\bullet`$

`"xo"` : number of crossovers

`"dxo"` : number of double crossovers

`"miss"` : number of missing values

The two statistics `"xo"` and `"dxo"` are obviously only useful for
constructed linkage maps. However, from my experience, they represent
the most vital two statistics for determining a linkage maps quality.
Inflated crossover or double crossover rates of any genotypes indicate a
lack of adherence to Mendelian genetics. For example, it can be quickly
calculated for a Doubled Haploid wheat population of size 100, one
recombination is approximately 1cM. Consequently, under Mendelian
genetics, a chromosome of length 200 cM would expect to have 200 random
crossovers with each genotype having an expected recombination rate of
2. Wheat is a hexaploid, so for excellent coverage over the whole genome
the expected recombination rate of any genotype is $`\sim`$42. This
number will obviously be reduced if the coverage of the genome is
incomplete. Significant recombination rates can be checked by manually
inputting a median recombination rate in the argument `xo.lambda` of
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md).
As an example, the genotype profiles of the `mapDH` can be displayed
using the command below and the resultant plot is given in the figure
below.

``` r

profileGen(mapDH, bychr = FALSE, stat.type = c("xo", "dxo", "miss"), id = "Genotype", xo.lambda = 25, layout = c(1,3), lty = 2)
```

![Genotype profiles of missing values, double recombinations and
recombinations for
\`mapDH\`.](diagnostics_files/figure-html/prof2-1.png)

Genotype profiles of missing values, double recombinations and
recombinations for `mapDH`.

The plot neatly displays, the number of missing values, double
crossovers and crossovers for each genotype in the order represented in
`mapDH`. The plot was aesthetically enhanced by including graphical
arguments `layout = c(1,3)` and `lty = 2` which are passed to the high
level lattice function
[`xyplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html). Plotting
these statistics simultaneously allows the users to quickly recognize
genotypes that may be problematic and worth further investigation. For
example, in the figure below the line DH218 was identified as having an
inflated recombination rate and removing it may improve the quality of
the linkage map.

The [`statMark()`](https://drj001.github.io/ASMap/reference/statMark.md)
and
[`profileMark()`](https://drj001.github.io/ASMap/reference/profileMark.md)
functions are commands for the calculation and profiling of statistics
associated with the markers or intervals of the linkage map.
[`profileMark()`](https://drj001.github.io/ASMap/reference/profileMark.md)
calls
[`statMark()`](https://drj001.github.io/ASMap/reference/statMark.md) and
graphically profiles user given statistics and will also return them
invisibly after plotting. The current marker statistics that can be
profiled are

$`\bullet`$

`"seg.dist"`: -log10 p-value from a test of segregation distortion

`"miss"`: proportion of missing values

`"prop"`: allele proportions

`"dxo"`: number of double crossovers

The current interval statistics that can be profiled are

$`\bullet`$

`"erf"`: estimated recombination fractions

`"lod"`: LOD score for the test of no linkage

`"dist"`: interval map distance

`"mrf"`: map recombination fraction

`"recomb"`: number of recombinations.

The function allows any or all of the statistics to be plotted
simultaneously on a multi-panel lattice display. This includes
combinations of marker and interval statistics. There is a `chr`
argument to subset the linkage map to user defined linkage groups. If
`crit.val = "bonf"` then markers that have significant segregation
distortion greater than the family wide alpha level of
0.05/no.of.markers will be annotated in marker panels. Similarly,
intervals that have a significantly weak linkage from a test of the
recombination fraction of $`r = 0.5`$ will also be annotated in the
interval panels. Similar to
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md),
graphical arguments can be used in
[`profileMark()`](https://drj001.github.io/ASMap/reference/profileMark.md)
and are passed to the high level lattice function
[`xyplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html). A plot of the
segregation distortion, double crossovers, estimated recombination
fractions and LOD scores for the whole genome of `mapDH` is given in the
figure below and created with the command

``` r

profileMark(mapDH, stat.type = c("seg.dist", "dxo", "erf", "lod"), id = "Genotype", layout = c(1,4), type = "l")
```

![Marker and interval profiles of segregation distortion, double
crossovers, estimated recombination fractions and LOD scores for
\`mapDH\`.](diagnostics_files/figure-html/prof3-1.png)

Marker and interval profiles of segregation distortion, double
crossovers, estimated recombination fractions and LOD scores for
`mapDH`.

All linkage groups are highlighted in a different colours to ensure they
can be identified clearly. The lattice panels ensure that marker and
interval statistics are seamlessly plotted together so problematic
regions or markers can be identified efficiently. In this command, the
`layout = c(1,4)` and `type = "l"` arguments are passed directly to the
high level [`xyplot()`](https://rdrr.io/pkg/lattice/man/xyplot.html)
lattice function to ensure a more aesthetically pleasing graphic.

### Genetic clones

In my experience, the assumptions of how the individuals of the
population are genetically related is rarely checked throughout the
construction process. Too often unconstructed or constructed linkage
maps contain individuals that are far too closely related beyond the
simple assumptions of the population. For example, in a DH population,
the assumption of independence would indicate that any two individuals
will, by chance, share half of their alleles. Any pairs of individuals
that significantly breach this assumption should be deemed suspicious
and queried.

The R/ASMap package contains a function for the detection and reporting
of the relatedness between individuals as well as a function for forming
consensus genotypes if genuine clones are found.

``` r

genClones(object, chr, tol = 0.9, id = "Genotype")
fixClones(object, gc, id = "Genotype", consensus = TRUE)
```

The
[`genClones()`](https://drj001.github.io/ASMap/reference/genClones.md)
function uses the power of
[`comparegeno()`](https://rdrr.io/pkg/qtl/man/comparegeno.html) from the
R/qtl package to perform the relatedness calculations. It then provides
a numerical breakdown of the relatedness between pairs of individuals
that share a proportion of alleles greater than `tol`. This breakdown
also includes the clonal group the pairs of individuals belong to. The
complete set of statistics allows users to make an informed decision
about the connectedness of the pairs of individuals. For example, using
the constructed map `mapDH` the list of clones can be found using

``` r

gc <- genClones(mapDH, tol = 0.9)
gc$cgd
```

         G1    G2 coef match diff na.both na.one group
    1  DH40  DH24    1   597    0       0      2     1
    2  DH59  DH53    1   597    0       0      2     2
    3  DH65  DH60    1   598    0       0      1     3
    4 DH143 DH139    1   597    0       0      2     4
    5 DH186 DH169    1   595    0       0      4     5

The reported information shows five pairs of clones that are, with the
exception of missing values, identical.

If clones are found then
[`fixClones()`](https://drj001.github.io/ASMap/reference/fixClones.md)
can be used to form consensus genotypes for each of the clonal groups.
By default it will intelligently collapse the allelic information of the
clones within each group (see the table below) to obtain a single
consensus genotype. Setting `consensus = FALSE` will choose a genotype
with the smallest proportion of missing values as the representative
genotype for the clonal group. In both cases genotypes names are given
an elongated name containing all genotype names of the clonal group
separated by an underscore.

``` r

mapDHg <- fixClones(mapDH, gc$cgd, consensus = TRUE)
levels(mapDHg$pheno[[1]])[grep("_", levels(mapDHg$pheno[[1]]))]
```

    [1] "DH139_DH143" "DH169_DH186" "DH24_DH40"   "DH53_DH59"   "DH60_DH65"  

| Genotype | M1  | M2  | M3  | M4  | M5  | M6  | M7  | M8  |
|:---------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| G1       | AA  | BB  | AA  | BB  | AA  | BB  | NA  | AA  |
| G2       | AA  | BB  | NA  | NA  | NA  | NA  | NA  | BB  |
| G3       | AA  | BB  | AA  | BB  | NA  | NA  | NA  | AA  |
| Cons.    | AA  | BB  | AA  | BB  | AA  | BB  | NA  | NA  |

Consensus genotype outcomes for 3 clones across 8 markers in a DH
population. {.table}
