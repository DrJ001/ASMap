# Package index

## Construction

Linkage group clustering, optimal marker ordering and genetic distance
estimation in a single pass, using the MSTmap algorithm. Both functions
are methods of the
[`mstmap()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
generic and share the full set of MSTmap parameters, differing only in
their input.

- [`mstmap(`*`<cross>`*`)`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
  :

  Extremely fast linkage map construction for qtl objects using MSTmap.

- [`mstmap(`*`<data.frame>`*`)`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
  : Extremely fast linkage map construction for data frame objects using
  MSTmap.

## Pre-construction

Diagnosis and curation of genotypes and markers before a map is built.

### Genotype diagnostics

- [`statGen()`](https://drj001.github.io/ASMap/reference/statGen.md) :
  Individual genotype statistics for an R/qtl cross object
- [`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md)
  : Profile individual genotype statistics for an R/qtl cross object
- [`genClones()`](https://drj001.github.io/ASMap/reference/genClones.md)
  : Find and report genotype clones
- [`fixClones()`](https://drj001.github.io/ASMap/reference/fixClones.md)
  : Consensus genotypes for clonal genotype groups

### Marker and interval diagnostics

- [`statMark()`](https://drj001.github.io/ASMap/reference/statMark.md) :
  Individual marker and interval statistics for an R/qtl cross object
- [`profileMark()`](https://drj001.github.io/ASMap/reference/profileMark.md)
  : Profile individual marker and interval statistics for an R/qtl cross
  object
- [`pValue()`](https://drj001.github.io/ASMap/reference/pValue.md) :
  P-value graph

### Setting markers aside

Pull markers of a given type out of the map and hold them in the cross
object, rather than discarding them permanently.

- [`pullCross()`](https://drj001.github.io/ASMap/reference/pullCross.md)
  : Pull markers from a linkage map.
- [`pp.init()`](https://drj001.github.io/ASMap/reference/pp.init.md) :
  Parameter initialization function

## Post-construction

Interpretation, display and manipulation of a constructed linkage map.

### Display

- [`heatMap()`](https://drj001.github.io/ASMap/reference/heatMap.md) :
  Heat map of the estimated pairwise recombination fractions and LOD
  linkage between markers.
- [`alignCross()`](https://drj001.github.io/ASMap/reference/alignCross.md)
  : Graphical linkage group identity and alignment.

### Reintroducing markers

- [`pushCross()`](https://drj001.github.io/ASMap/reference/pushCross.md)
  : Push markers into an established R/qtl linkage map.

### Map manipulation

- [`breakCross()`](https://drj001.github.io/ASMap/reference/breakCross.md)
  :

  Break linkage groups of an qtl cross object

- [`mergeCross()`](https://drj001.github.io/ASMap/reference/mergeCross.md)
  :

  Merge linkage groups of an qtl cross object

- [`combineMap()`](https://drj001.github.io/ASMap/reference/combineMap.md)
  :

  Combine linkage maps from multiple qtl cross objects

- [`subsetCross()`](https://drj001.github.io/ASMap/reference/subsetCross.md)
  : Subset an R/qtl object

- [`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md) :
  Very quick estimation of genetic map distances.

## Datasets

Constructed and unconstructed marker sets used throughout the
documentation. These are not lazy-loaded, so an explicit
[`data()`](https://rdrr.io/r/utils/data.html) call is required.

- [`mapDH`](https://drj001.github.io/ASMap/reference/mapDH.md) : A
  constructed linkage map for a doubled haploid wheat population
- [`mapDHf`](https://drj001.github.io/ASMap/reference/mapDHf.md) : An
  unconstructed marker set for a doubled haploid wheat population
- [`mapBC`](https://drj001.github.io/ASMap/reference/mapBC.md) : A
  constructed linkage map for a backcross barley population
- [`mapBCu`](https://drj001.github.io/ASMap/reference/mapBCu.md) : An
  unconstructed marker set for a backcross barley population
- [`mapF2`](https://drj001.github.io/ASMap/reference/mapF2.md) :
  Simulated constructed linkage map for a self pollinated F2 barley
  population

## Package

- [`ASMap-package`](https://drj001.github.io/ASMap/reference/ASMap-package.md)
  : Additional functions for linkage map construction and manipulation
  of R/qtl objects.
