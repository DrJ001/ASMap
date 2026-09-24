# ASMap - Accurate and Speedy Linkage Map Construction

<!-- badges: start -->
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.0.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN status](https://www.r-pkg.org/badges/version/ASMap)](https://CRAN.R-project.org/package=ASMap)
[![R-CMD-check](https://github.com/DrJ001/ASMap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DrJ001/ASMap/actions/workflows/R-CMD-check.yaml)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/ASMap)](https://CRAN.R-project.org/package=ASMap)
<!-- badges: end -->

**Authors**: Julian Taylor & David Butler

---

## Overview

**ASMap** is a linkage map construction, manipulation and diagnosis toolkit for
biparental populations, built around the efficient **MSTmap** algorithm of
[Wu *et al.* (2008)](https://doi.org/10.1371/journal.pgen.1000212). MSTmap uses the
minimum spanning tree of a graph to cluster markers into linkage groups **and**
optimally order the markers within each group, in a single pass.

| Pillar | Functions | Description |
|---|---|---|
| Construction | `mstmap()` | Simultaneous linkage group clustering, optimal marker ordering and genetic distance estimation, for `data.frame` or R/qtl `cross` input |
| Marker pull/push | `pullCross()` • `pushCross()` | Temporarily set problematic markers aside before construction and reintroduce them afterwards, rather than discarding them |
| Diagnosis | `profileGen()` • `profileMark()` • `heatMap()` • `alignCross()` | Numerical and graphical interrogation of genotypes, markers and intervals, before and after construction |

Objects use the **R/qtl** `cross` format throughout, so ASMap and
[qtl](https://CRAN.R-project.org/package=qtl) functions can be used synergistically
on the same object at any stage.

Map construction with ASMap follows the same three-stage workflow: **pre-construction**
diagnosis and marker curation, **construction**, and **post-construction**
diagnosis and refinement. The three parts below mirror that sequence.

### Supported populations

| `mstmap` `pop.type` | R/qtl class | Population |
|---|---|---|
| `"BC"` | `"bc"` | Backcross |
| `"DH"` | `"dh"` | Doubled haploid |
| `"ARIL"` | `"riself"` | Advanced recombinant inbred (treated as biparental) |
| `"RILn"` | `"bcsft"` | Non-advanced RIL with *n* generations of selfing (selfed F2, F3, …, F*n*) |

Doubled haploid populations read in with `read.cross()` inherit class `"bc"`; assigning
class `"dh"` gives equivalent construction results. `"RILn"` populations are imported
as `"f2"` and must be converted with `qtl::convert2bcsft(F.gen = n, BC.gen = 0)`.

### Performance

Full construction from scratch - clustering, ordering and distance estimation
(`bychr = FALSE`):

| Dataset | Individuals | Markers | Time |
|---|---|---|---|
| `mapDH` | 218 | 599 | 0.5 s |
| `mapF2` | 250 | 700 | 6.3 s |
| `mapBC` | 300 | 3,019 | 17.3 s |

<sub>Single run, R 4.4.3 on 64-bit Windows. Indicative only - timings depend on the
machine, population type and marker data.</sub>

---

## Requirements

ASMap requires **R >= 3.0.0** and depends on
[qtl](https://CRAN.R-project.org/package=qtl) and
[lattice](https://CRAN.R-project.org/package=lattice), importing
[fields](https://CRAN.R-project.org/package=fields),
[RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer) and
[gtools](https://CRAN.R-project.org/package=gtools). All are available from CRAN.

No external software or licence is needed. The MSTmap C++ source code is bundled with
the package and compiled on installation.

---

## Installation

Install the released version from CRAN:

```r
install.packages("ASMap")
```

Or the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("DrJ001/ASMap")
```

---

## Articles

Full documentation is published at
**[drj001.github.io/ASMap](https://drj001.github.io/ASMap/)**, as a set of
articles covering every function in the package, the technical basis of the MSTmap
algorithm, and a worked example taking an unconstructed 3,023 marker barley backcross
through to a diagnosed seven linkage group map.

#### Getting started

| Article | Description |
|---|---|
| [Getting started with ASMap](https://drj001.github.io/ASMap/articles/getting-started.html) | Construct a linkage map from a set of marker scores, the choice of `p.value`, and an orientation to the remaining documentation. |

#### Function guides

| Article | Description |
|---|---|
| [Constructing a linkage map](https://drj001.github.io/ASMap/articles/constructing-a-map.html) | The two construction functions, the population types and allele codings they accept, and the four modes of reconstruction. |
| [Pulling and pushing markers](https://drj001.github.io/ASMap/articles/pulling-and-pushing.html) | Setting problematic markers aside before construction and reintroducing them afterwards, and the reversal of the thresholds between the two. |
| [Diagnosing genotypes and markers](https://drj001.github.io/ASMap/articles/diagnostics.html) | Numerical and graphical profiling of genotypes, markers and intervals, and the detection of genetic clones. |
| [Heat maps](https://drj001.github.io/ASMap/articles/heat-maps.html) | Displaying pairwise recombination fractions and LOD scores together, and matching the two scales. |
| [Manipulating linkage maps](https://drj001.github.io/ASMap/articles/manipulating-maps.html) | Breaking, merging, combining and subsetting maps, and rapid estimation of genetic distances. |

#### Worked example

| Article | Description |
|---|---|
| [Worked example I: construction](https://drj001.github.io/ASMap/articles/worked-example-construction.html) | A barley backcross taken from 3,023 unassigned markers through pre-construction diagnosis to a constructed map of nine linkage groups. |
| [Worked example II: refinement](https://drj001.github.io/ASMap/articles/worked-example-refinement.html) | Reintroducing the markers held aside, merging to the expected seven linkage groups, and post-construction development including fine mapping. |

#### Background

| Article | Description |
|---|---|
| [How the MSTmap algorithm works](https://drj001.github.io/ASMap/articles/mstmap-algorithm.html) | The minimum spanning tree formulation, clustering, marker ordering, imputation of missing scores and detection of genotyping errors. |
| [Notes on the MSTmap algorithm](https://drj001.github.io/ASMap/articles/algorithm-notes.html) | How genetic distances become inflated, and what `mvest.bc` and `detectBadData` do in practice. |

The documentation is summarised in the *Journal of Statistical Software* paper
[Taylor & Butler (2017)](https://doi.org/10.18637/jss.v079.i06).

---

## Part 1 - Pre-construction: Diagnosis and Marker Curation

A linkage map should not be built by throwing raw marker data at a construction
algorithm. The functions below support a pre-construction checklist: check missing
values across genotypes and markers, check for genetic clones, check for excessive
segregation distortion, and set aside markers that are co-located or otherwise
problematic.

### Genotype statistics - `statGen()` and `profileGen()`

`statGen()` returns per-genotype statistics; `profileGen()` plots the same statistics
as a multi-panel lattice display and returns them invisibly.

```r
statGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"),
        id = "Genotype")

profileGen(cross, chr, bychr = TRUE, stat.type = c("xo", "dxo", "miss"),
           id = "Genotype", xo.lambda = NULL, ...)
```

| Argument | Description |
|---|---|
| `cross` | An R/qtl `cross` object inheriting `"bc"`, `"dh"`, `"riself"` or `"bcsft"` |
| `chr` | Character vector of linkage group names used to subset the map |
| `bychr` | If `TRUE` (default), statistics are returned/plotted by linkage group; if `FALSE`, across the whole genome |
| `stat.type` | Any combination of `"xo"` (number of crossovers), `"dxo"` (double crossovers) and `"miss"` (missing values) |
| `id` | Column of `cross$pheno` containing the genotype names. Default `"Genotype"` |
| `xo.lambda` | `profileGen()` only. Expected genome-wide recombination rate; genotypes whose crossover count significantly exceeds it are annotated on the plot |
| `...` | `profileGen()` only. Passed to the high level lattice plot |

When `xo.lambda` is supplied, each genotype's recombination count is tested against
that expected rate with a one-tailed Poisson test, and lines significant at a
Bonferroni-adjusted 0.05 level are flagged. The returned object gains a logical
`xo.lambda` element, which pairs directly with `subsetCross()` to drop the offenders.

### Genetic clones - `genClones()` and `fixClones()`

Highly related or duplicated individuals inflate segregation distortion, so they
should be resolved before construction.

```r
genClones(object, chr, tol = 0.9, id = "Genotype")

fixClones(object, gc, id = "Genotype", consensus = TRUE)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl cross object of any class structure |
| `chr` | Character vector of linkage group names |
| `tol` | Pairs of genotypes sharing a proportion of matching alleles above this value are reported. Default `0.9` |
| `id` | Column of `object$pheno` containing the genotype names. Default `"Genotype"` |
| `gc` | `fixClones()` only. The data frame of clone information, normally `genClones()$cgd` |
| `consensus` | `fixClones()` only. If `TRUE` (default), a consensus genotype is formed for each clonal group by intelligently collapsing alleles; if `FALSE`, only the group member with the fewest missing alleles is retained |

`genClones()` extends `qtl::comparegeno()` by returning breakdown statistics for each
reported pair in `$cgd`, including a `group` column identifying the clonal group.
**Inspect `$cgd` before acting on it** - pairs with few co-scored markers are weak
evidence of cloning and are best removed from the table before calling `fixClones()`.
`fixClones()` can be used at any stage, as it retains linkage group and marker
position information.

### Marker and interval statistics - `statMark()` and `profileMark()`

```r
statMark(cross, chr, stat.type = c("marker", "interval"),
         map.function = "kosambi")

profileMark(cross, chr, stat.type = "marker", use.dist = TRUE,
            map.function = "kosambi", crit.val = NULL,
            display.markers = FALSE, mark.line = FALSE, ...)
```

| Argument | Description |
|---|---|
| `cross` | An R/qtl `cross` object inheriting `"bc"`, `"dh"`, `"riself"` or `"bcsft"` |
| `chr` | Character vector of linkage group names used to subset the map |
| `stat.type` | `"marker"` and/or `"interval"`. `profileMark()` also accepts individual statistics: marker - `"seg.dist"`, `"miss"`, `"prop"`, `"dxo"`; interval - `"erf"`, `"lod"`, `"dist"`, `"mrf"`, `"recomb"` |
| `map.function` | `"kosambi"` (default), `"haldane"`, `"morgan"` or `"cf"`, used for interval statistics |
| `use.dist` | `profileMark()` only. If `TRUE` (default) markers are positioned by actual map distance, otherwise equidistantly |
| `crit.val` | `profileMark()` only. Set to `"bonf"` to annotate markers whose segregation distortion p-value, and intervals whose test of no linkage, fall below a Bonferroni-adjusted 0.05 level |
| `display.markers` | `profileMark()` only. Display marker names on the bottom axis. Default `FALSE` |
| `mark.line` | `profileMark()` only. Draw vertical lines at marker positions. Default `FALSE` |
| `...` | `profileMark()` only. Passed to the high level lattice plot |

`"marker"` statistics use `qtl::geno.table()` to give segregation distortion, allele
and missing value proportions, and double crossover counts at each marker under the
current order. `"interval"` statistics give the estimated recombination fraction and
LOD score between adjacent markers, map distances, and actual recombination counts.

### Choosing `p.value` - `pValue()`

The `p.value` argument of `mstmap()` governs how markers are split into linkage
groups, and the appropriate value depends strongly on population size - larger
populations require a smaller `p.value`. `pValue()` plots the relationship so an
informed choice can be made.

```r
pValue(dist = seq(25, 40, by = 5), pop.size = 100:500,
       map.function = "kosambi", LOD = FALSE)
```

| Argument | Description |
|---|---|
| `dist` | Numeric range of genetic distances in cM |
| `pop.size` | Numeric range of population sizes |
| `map.function` | `"kosambi"` (default), `"haldane"`, `"morgan"` or `"cf"` |
| `LOD` | If `TRUE` the LOD score of linkage is calculated; if `FALSE` (default) the minus log10 p-value used to threshold the Hoeffding inequality |

### Setting markers aside - `pullCross()` and `pp.init()`

Pruning before construction is common, and usually permanent. `pullCross()` instead
*removes markers to a holding element* of the cross object, retaining everything
needed to push them back later.

```r
pullCross(object, chr, type = c("co.located", "seg.distortion", "missing"),
          pars = NULL, replace = FALSE, ...)

pp.init(seg.thresh = 0.05, seg.ratio = NULL, miss.thresh = 0.1,
        max.rf = 0.25, min.lod = 3)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl `cross` object inheriting `"bc"`, `"dh"`, `"riself"` or `"bcsft"` |
| `chr` | Character vector of linkage group names to subset the map before pulling |
| `type` | One of `"co.located"`, `"seg.distortion"` or `"missing"` |
| `pars` | List of threshold parameters; `NULL` (default) calls `pp.init()` with its defaults |
| `replace` | If `TRUE`, pulled markers replace what already resides in the `type` element. Default `FALSE` |
| `...` | Currently ignored |

Thresholds supplied through `pars`:

| Parameter | Role in `pullCross()` | Default |
|---|---|---|
| `seg.thresh` | Markers are pulled if their segregation distortion p-value is **less** than this | `0.05` |
| `seg.ratio` | Alternative to `seg.thresh`; character ratio `"AA:BB"` or `"AA:AB:BB"` | `NULL` |
| `miss.thresh` | Markers are pulled if their missing value proportion is **greater** than this | `0.1` |
| `max.rf` | Maximum recombination fraction used when clustering markers back in (`pushCross()`) | `0.25` |
| `min.lod` | Minimum LOD score used when clustering markers back in (`pushCross()`) | `3` |

For `type = "co.located"` the first marker of each co-locating group is retained in
the map as a reference marker and the rest are set aside, keyed to it.

---

## Part 2 - Construction functions

Both construction functions are methods of the generic `mstmap()` and share the
complete set of MSTmap algorithm parameters. They differ in their input: a raw data
frame of marker scores, or an existing R/qtl `cross` object.

### From a data frame - `mstmap.data.frame()`

Use this when marker data has not yet been imported into an R/qtl object. The format
follows the marker file required by the standalone MSTmap software: **markers in rows,
genotypes in columns**, marker names in `rownames`, genotype names in `names`, and
every column of class `"character"` (not factors).

```r
mstmap(object, pop.type = "DH", dist.fun = "kosambi", objective.fun = "COUNT",
       p.value = 1e-06, noMap.dist = 15, noMap.size = 0, miss.thresh = 1,
       mvest.bc = FALSE, detectBadData = FALSE, as.cross = TRUE,
       return.imputed = FALSE, trace = FALSE, ...)
```

| Argument | Description |
|---|---|
| `object` | A `data.frame` of marker scores, markers in rows and genotypes in columns |
| `pop.type` | `"DH"` (default), `"BC"`, `"ARIL"` or `"RILn"` |
| `as.cross` | If `TRUE` (default) an R/qtl cross object of the appropriate class is returned; if `FALSE` a `data.frame` with linkage group, marker position and distance columns |

Allele coding is strict. For `"BC"`, `"DH"` and `"ARIL"` the two alleles must be
`"A"`/`"a"` and `"B"`/`"b"`; `"ARIL"` additionally assumes heterozygotes have been set
to missing. For `"RILn"`, phase-unknown heterozygotes must be coded `"X"`. Missing
scores are `"U"` or `"-"` for all population types.

### From an R/qtl cross object - `mstmap.cross()`

The more flexible of the two, and the recommended entry point. It accepts a
constructed *or* unconstructed cross object and can rebuild the whole map or re-order
markers within nominated linkage groups.

```r
mstmap(object, chr, id = "Genotype", bychr = TRUE, suffix = "numeric",
       anchor = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
       p.value = 1e-06, noMap.dist = 15, noMap.size = 0, miss.thresh = 1,
       mvest.bc = FALSE, detectBadData = FALSE, return.imputed = FALSE,
       trace = FALSE, ...)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl `cross` object inheriting `"bc"`, `"dh"`, `"riself"` or `"bcsft"` |
| `chr` | Character vector of linkage group names requiring reconstruction and/or re-ordering |
| `id` | Column of `object$pheno` uniquely identifying genotype names. Default `"Genotype"` |
| `bychr` | If `TRUE` (default), split linkage groups if required and order markers within them; if `FALSE`, bulk all linkage groups and reconstruct from scratch |
| `suffix` | `"numeric"` or `"alpha"` - whether numeric or alphabetic values are appended to linkage group names when groups are split |
| `anchor` | If `TRUE`, the inputted marker order is respected regardless of `chr` and `bychr`. Default `FALSE` |

Parameters shared by both methods:

| Argument | Description |
|---|---|
| `dist.fun` | Distance function for genetic distances: `"kosambi"` (default) or `"haldane"` |
| `objective.fun` | `"COUNT"` (default) minimises the sum of recombination events between markers; `"ML"` maximises the likelihood objective function |
| `p.value` | Threshold for clustering markers into linkage groups. Default `1e-06`. **A value greater than 1 turns clustering off**, and all inputted markers are assumed to belong to the same linkage group |
| `noMap.dist` | Smallest genetic distance at which a set of isolated markers appears distinct from other linked markers. Default `15` |
| `noMap.size` | Maximum size of the isolated marker linkage groups identified by `noMap.dist`. Set to `0` (default) to turn this off |
| `miss.thresh` | Markers with a proportion of missing scores above this are excluded from the map. Default `1` |
| `mvest.bc` | If `TRUE`, missing markers are imputed before clustering. Restricted to BC/DH/ARIL populations. Default `FALSE` |
| `detectBadData` | If `TRUE`, likely genotyping errors are detected, set to missing and imputed during ordering, and reported in the `trace` file. Restricted to BC/DH/ARIL populations. Default `FALSE` |
| `return.imputed` | If `TRUE`, the imputed marker probability matrix is returned for the constructed linkage groups. Default `FALSE` |
| `trace` | If `FALSE` (default) minimal output is piped to the screen; if `TRUE`, detailed MSTmap output is written to `"MSToutput.txt"`, equivalent to running the MSTmap executable from the command line |
| `...` | Currently ignored |

Two settings of `p.value` are worth knowing as idioms:

```r
# Cluster and order from scratch, ignoring existing linkage groups
map <- mstmap(cross, bychr = FALSE, dist.fun = "kosambi", p.value = 1e-12)

# Re-order within existing linkage groups only, never splitting them
map <- mstmap(cross, bychr = TRUE, dist.fun = "kosambi", p.value = 2)
```

---

## Part 3 - Post-construction: Interpretation, Display and Manipulation

### Heat maps - `heatMap()`

The single most informative check on a constructed map. `heatMap()` displays pairwise
estimated recombination fractions and pairwise LOD scores of linkage **with a separate
legend for each**, so the two scales can be tuned until their heat matches. Consistent
heat within linkage groups indicates strong linkage between nearby markers; blocks of
shared heat between groups indicate groups that should be merged.

```r
heatMap(x, chr, mark, what = c("both", "lod", "rf"), lmax = 12, rmin = 0,
        markDiagonal = FALSE,
        color = rev(colorRampPalette(brewer.pal(11, "Spectral"))(256)), ...)
```

| Argument | Description |
|---|---|
| `x` | An R/qtl `cross` object |
| `chr` | Character vector of linkage group names to subset the object |
| `mark` | Subsets linkage groups further into marker subsets; a numeric vector applied to all groups, or a list of numeric vectors named by linkage group |
| `what` | `"both"` (default) plots LOD on the lower triangle and recombination fractions on the upper; `"lod"` or `"rf"` plot one only |
| `lmax` | LOD scores above this threshold are plotted in the same colour. Default `12` |
| `rmin` | Recombination fractions below this threshold are plotted in the same colour. Default `0` |
| `markDiagonal` | If `TRUE`, borders are drawn around the diagonal elements. Default `FALSE` |
| `color` | Colour spectrum, defaulting to a reversed `RColorBrewer` `"Spectral"` palette |
| `...` | Additional arguments, largely passed to `fields::image.plot()`. **Avoid `legend.args`**, which is used internally |

### Linkage group identity and orientation - `alignCross()`

Newly constructed linkage groups carry arbitrary names and orientation. `alignCross()`
plots marker positions against one or more reference maps to identify which group is
which and whether it needs inverting - a positive slope indicates correct
orientation, a negative slope indicates the group should be reversed.

```r
alignCross(object, chr, maps, ...)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl cross object of any class structure |
| `chr` | Character vector of linkage group names, or a logical vector of length equal to the number of linkage groups |
| `maps` | A **named list** of cross objects or data frames containing markers present in `object` |
| `...` | Passed to the high level lattice plot |

Data frame elements of `maps` must contain the columns `"marker"`, `"ref.chr"` and
`"ref.dist"`, or an error is produced.

### Reintroducing markers - `pushCross()`

The complement of `pullCross()`. Markers held aside are pushed back into the
established map, then the map is re-ordered by linkage group.

```r
pushCross(object, type = c("co.located", "seg.distortion", "missing", "unlinked"),
          unlinked.chr = NULL, pars = NULL, ...)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl `cross` object inheriting `"bc"`, `"dh"`, `"riself"` or `"bcsft"` |
| `type` | One of `"co.located"`, `"seg.distortion"`, `"missing"` or `"unlinked"` |
| `unlinked.chr` | Character vector of linkage group names holding markers to push into the remaining groups. Used only with `type = "unlinked"` |
| `pars` | List of threshold parameters; `NULL` (default) calls `pp.init()` with its defaults |
| `...` | Currently ignored |

The threshold logic is **deliberately the reverse** of `pullCross()`:

| Parameter | Role in `pushCross()` | Default |
|---|---|---|
| `seg.thresh` | Markers are pushed back if their distortion p-value is **greater** than this | `0.05` |
| `miss.thresh` | Markers are pushed back if their missing value proportion is **less** than this | `0.1` |

So markers are typically pulled at one threshold and only pushed back if they clear a
more permissive one, letting a marker earn its place in the map.

`type = "unlinked"` is the mechanism for pushing new markers into an established map -
the basis of the fine mapping workflow in the vignette. For all types other than
`"co.located"`, markers are allocated to linkage groups by a fast clustering method
that reduces the map to a skeleton marker set, tuned via `max.rf` and `min.lod`.
Avoid `"UL"` in linkage group names, as it is used internally for unlinked groups.

> `pushCross()` allocates markers to linkage groups but does **not** re-construct the
> map. Follow it with `mstmap(..., bychr = TRUE, p.value = 2)` to place the new
> markers optimally.

### Map manipulation - `breakCross()`, `mergeCross()`, `combineMap()`, `subsetCross()`

```r
breakCross(cross, split = NULL, suffix = "numeric", sep = ".")
mergeCross(cross, merge = NULL, gap = 5)
combineMap(..., id = "Genotype", keep.all = TRUE, merge.by = "genotype")
subsetCross(cross, chr, ind, ...)
```

| Function | Key arguments |
|---|---|
| `breakCross()` | `split` - a list named by linkage group containing the marker names *immediately preceding* each split point. `suffix` - `"numeric"` or `"alpha"`, or a list of explicit new names. `sep` - separator between old name and suffix, default `"."` |
| `mergeCross()` | `merge` - a list whose elements contain the linkage groups to merge, each element named by the proposed new linkage group name. `gap` - cM gap inserted between merged elements, default `5` |
| `combineMap()` | `...` - any number of cross objects, all of identical class. `id` - common genotype column, default `"Genotype"`. `keep.all` - retain genotypes absent from some maps, default `TRUE`. `merge.by` - `"genotype"` (default) or `"marker"` |
| `subsetCross()` | `chr` - logical, numeric or character vector of linkage groups. `ind` - logical or numeric vector of individuals |

Two cautions worth repeating:

- **Use `subsetCross()`, not `subset.cross()`**, whenever the object carries
  `"co.located"`, `"seg.distortion"` or `"missing"` elements from a `pullCross()` call.
  `subsetCross()` subsets those elements too and recalculates their statistics;
  `subset.cross()` silently leaves them inconsistent with the map.
- `combineMap()` applies **no consensus map algorithm**. Chromosome identity and
  genetic distances are taken from the first appearance of a marker across the maps
  supplied, and the result is not re-constructed - run `mstmap()` afterwards.

### Fast distance estimation - `quickEst()`

```r
quickEst(object, chr, map.function = "kosambi", ...)
```

| Argument | Description |
|---|---|
| `object` | An R/qtl cross object of any class structure |
| `chr` | Character vector of linkage group names requiring (re)estimation of distances |
| `map.function` | `"kosambi"` (default), `"haldane"`, `"morgan"` or `"cf"` |
| `...` | Passed to `qtl::argmax.geno()` |

`quickEst()` sidesteps the slow hidden Markov algorithm of `qtl::est.map()`. Initial
conservative distances are obtained by inverting recombination fractions from
`est.rf()`, then passed to `argmax.geno()`, where missing allele scores are imputed by
the Viterbi algorithm and distances re-estimated.

---

## A worked example

<details>
<summary><b>Click to expand: the condensed barley backcross workflow from the vignette</b></summary>

Starting from `mapBCu` - 3,023 markers scored on 326 individuals, with no linkage
groups assigned.

### 1. Pre-construction diagnosis

```r
library(ASMap)
data(mapBCu)

# Drop genotypes with excessive missing data
sg <- statGen(mapBCu, bychr = FALSE, stat.type = "miss")
mapBC1 <- subset(mapBCu, ind = sg$miss < 1600)

# Identify genetic clones, then form consensus genotypes. Inspect gc$cgd first -
# groups with few co-scored markers are not convincing evidence of cloning.
gc <- genClones(mapBC1, tol = 0.95)
mapBC2 <- fixClones(mapBC1, gc$cgd[-c(1, 4, 5), ], consensus = TRUE)

# Profile segregation distortion, allele proportions and missingness
profileMark(mapBC2, stat.type = c("seg.dist", "prop", "miss"),
            crit.val = "bonf", layout = c(1, 4), type = "l", cex = 0.5)

# Drop the handful of extremely distorted markers
mm <- statMark(mapBC2, stat.type = "marker")$marker$AB
mapBC3 <- drop.markers(mapBC2, markernames(mapBC2)[mm > 0.98 | mm < 0.2])
```

### 2. Set problematic markers aside

```r
mapBC3 <- pullCross(mapBC3, type = "missing", pars = list(miss.thresh = 0.1))
mapBC3 <- pullCross(mapBC3, type = "seg.distortion", pars = list(seg.thresh = "bonf"))
mapBC3 <- pullCross(mapBC3, type = "co.located")

sum(nmar(mapBC3))   # 847 markers set aside
#> [1] 2173
```

### 3. Construct

```r
mapBC4 <- mstmap(mapBC3, bychr = FALSE, dist.fun = "kosambi", p.value = 1e-12)

nchr(mapBC4)
#> [1] 9

heatMap(mapBC4, lmax = 70)
```

### 4. Remove problem genotypes and rebuild

Linkage groups beyond ~200 cM suggest genotypes with inflated recombination counts.

```r
pg <- profileGen(mapBC4, bychr = FALSE, stat.type = c("xo", "dxo", "miss"),
                 id = "Genotype", xo.lambda = 14, layout = c(1, 3), lty = 2, cex = 0.7)

# subsetCross() keeps the pulled-marker elements in step
mapBC5 <- subsetCross(mapBC4, ind = !pg$xo.lambda)
mapBC6 <- mstmap(mapBC5, bychr = TRUE, dist.fun = "kosambi", p.value = 1e-12)
```

### 5. Push markers back, merge and refine

```r
# Reintroduce the markers with 10-20% missing values
mapBC6 <- pushCross(mapBC6, type = "missing",
                    pars = list(miss.thresh = 0.22, max.rf = 0.3))

# The heat map now shows genuine linkage between L.3/L.5 and L.8/L.9, giving
# the 7 linkage groups expected for the barley genome
mapBC6 <- mergeCross(mapBC6, merge = list("L.3" = c("L.3", "L.5"),
                                          "L.8" = c("L.8", "L.9")))
names(mapBC6$geno) <- paste("L.", 1:7, sep = "")

mapBC7 <- mstmap(mapBC6, bychr = TRUE, dist.fun = "kosambi", p.value = 2)
chrlen(mapBC7)
```

The vignette continues from here into fine mapping, combining maps of the same
population, and resolving unknown linkage groups.

</details>

---

## Datasets

| Dataset | Class | Population | Individuals | Markers | Description |
|---|---|---|---|---|---|
| `mapDH` | `bc` | Doubled haploid (wheat) | 218 | 599 | Constructed map, 23 linkage groups |
| `mapDHf` | `data.frame` | Doubled haploid (wheat) | 218 | 599 | The same marker set, unconstructed, formatted for `mstmap.data.frame()` |
| `mapBC` | `bc` | Backcross (barley) | 300 | 3,019 | Constructed map, 7 linkage groups |
| `mapBCu` | `bc` | Backcross (barley) | 326 | 3,023 | Unconstructed marker set, all markers in one group |
| `mapF2` | `bcsft` | Selfed F2 (barley) | 250 | 700 | Simulated constructed map, 7 linkage groups |

These are not lazy-loaded, so an explicit `data()` call is required:

```r
data(mapDH)
```

---

## Citation

If you find this package useful, please cite it. Run `citation("ASMap")`, or:

> Taylor, J. & Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage Map
> Construction and Diagnosis. *Journal of Statistical Software*, **79**(6), 1-29.
> [doi:10.18637/jss.v079.i06](https://doi.org/10.18637/jss.v079.i06)

---

## References

Wu, Y., Bhat, P., Close, T.J. & Lonardi, S. (2008) Efficient and Accurate
Construction of Genetic Linkage Maps from the Minimum Spanning Tree of a Graph.
*PLoS Genetics*, **4**(10), e1000212.

Taylor, J. & Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage Map
Construction and Diagnosis. *Journal of Statistical Software*, **79**(6), 1-29.

Broman, K.W. & Sen, S. (2009) *A Guide to QTL Mapping with R/qtl*. Springer, New York.

---

## Getting help

Found a bug, or have a question about map construction?
[Open an issue](https://github.com/DrJ001/ASMap/issues).

## License

GPL (>= 2). The bundled MSTmap C++ source code is by Yonghui Wu, Timothy Close and
Stefano Lonardi.
