# Articles

### Getting started

- [Getting started with
  ASMap](https://drj001.github.io/ASMap/articles/getting-started.md):

  Construct a linkage map from a set of marker scores, and an
  orientation to the remainder of the documentation.

### Function guides

- [Constructing a linkage
  map](https://drj001.github.io/ASMap/articles/constructing-a-map.md):

  The two construction functions, for a data frame of marker scores and
  for an R/qtl cross object, and the MSTmap parameters they share.

- [Pulling and pushing
  markers](https://drj001.github.io/ASMap/articles/pulling-and-pushing.md):

  Setting problematic markers aside before construction and
  reintroducing them afterwards, in place of their permanent removal.

- [Diagnosing genotypes and
  markers](https://drj001.github.io/ASMap/articles/diagnostics.md):

  Numerical and graphical profiling of individual genotypes, markers and
  intervals, and the identification of genetic clones.

- [Heat maps](https://drj001.github.io/ASMap/articles/heat-maps.md):

  Displaying pairwise recombination fractions and LOD scores of linkage
  together, with a separate legend for each.

- [Manipulating linkage
  maps](https://drj001.github.io/ASMap/articles/manipulating-maps.md):

  Breaking, merging, subsetting and combining linkage maps, and fast
  estimation of genetic distances.

### Worked example

- [Worked example I:
  construction](https://drj001.github.io/ASMap/articles/worked-example-construction.md):

  A barley backcross marker set taken from 3,023 unassigned markers
  through pre-construction diagnosis to a constructed map of nine
  linkage groups.

- [Worked example II:
  refinement](https://drj001.github.io/ASMap/articles/worked-example-refinement.md):

  Continuing the barley backcross example: pushing the held-back markers
  into the map, merging linkage groups to the expected seven, and
  post-construction development including fine mapping.

### Background

- [How the MSTmap algorithm
  works](https://drj001.github.io/ASMap/articles/mstmap-algorithm.md):

  The minimum spanning tree formulation underlying ASMap: the clustering
  of markers into linkage groups, the determination of an optimal order
  within a group, and the imputation of missing scores and detection of
  genotyping errors that are integrated with it.

- [Notes on the MSTmap
  algorithm](https://drj001.github.io/ASMap/articles/algorithm-notes.md):

  How genetic distances are calculated, and what the mvest.bc and
  detectBadData arguments actually do.
