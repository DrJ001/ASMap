# Additional functions for linkage map construction and manipulation of R/qtl objects.

Additional functions for linkage map construction and manipulation of
R/qtl objects. This includes extremely fast linkage map clustering and
marker ordering using MSTmap (see Wu et al., 2008).

## Details

|          |            |
|----------|------------|
| Package: | ASMap      |
| Type:    | Package    |
| Version: | 1.0-4      |
| Date:    | 2018-10-24 |
| License: | GPL 2      |

Welcome to the ASMap package!

One of the fundamental reasons why this package exists was to utilize
and implement the source code for the the Minimum Spanning Tree
algorithm derived in Wu et al. (2008) (reference below) for linkage map
construction. The algorithm is lightning quick at linkage group
clustering and optimal marker ordering and can handle large numbers of
markers.

The package contains two very efficient functions, `mstmap.data.frame`
and `mstmap.cross`, that provide users with a highly flexible set
linkage map construction methods using the MSTmap algorithm.
`mstmap.data.frame` constructs a linkage map from a data frame of
genetic marker data and will use the entire contents of the object to
form linkage groups and optimally order markers within each linkage
group. `mstmap.cross` is a linkage map construction function for qtl
package objects and can be used to construct linkage maps in a flexible
number of ways. See
[`?mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
for complete details.

To complement the computationally efficient linkage map construction
functions, the package also contains functions `pullCross` and
`pushCross` that allow the pulling/pushing markers of different types to
and from the linkage map. This system gives users the ability to
initially pull markers aside that are not needed for immediate
construction and push them back later if required. There are also
functions for fast numerical and graphical diagnosis of unconstructed
and constructed linkage maps. Specifically, there is an improved
`heatMap` that graphically displays pairwise recombination fractions and
LOD scores with separate legends for each. `profileGen` can be used to
simultaneously profile multiple statistics such as recombination counts
and double recombination counts for individual lines across the
constructed linkage map. `profileMark` allows simultaneous graphical
visualization of marker or interval statistics profiles across the
genome or subsetted for a predefined set of linkage groups. Graphical
identification and orientation of linkage groups using reference linkage
maps can be conducted using `alignCross`. All of these graphical
functions utilize the power of the advanced graphics package lattice to
provide seamless multiple displays.

Other miscellaneous utilities for qtl objects include

- `mergeCross`: Merging of linkage groups

- `breakCross`: Breaking of linkage groups

- `combineMap`: Combining linkage maps

- `quickEst`: Very quick estimation of genetic map distances

- `genClones`: Reporting genotype clones

- `fixClones`: Consensus genotypes for clonal groups

Complete documentation for the package is published at
<https://drj001.github.io/ASMap/>, as a set of articles containing
detailed explanations of each function and how it may be used to perform
efficient map construction. These include a fully worked example
involving pre-construction diagnostics, linkage map construction and
post-construction diagnostics, and extending to post-construction
techniques such as fine mapping and the combining of linkage maps.

A short introduction is installed with the package and may be viewed
with
[`vignette("ASMap")`](https://drj001.github.io/ASMap/articles/ASMap.md).
The package is succinctly described in the Journal of Statistical
Software publication Taylor and Butler (2017) referenced below.

## Author

Julian Taylor, Dave Butler, Timothy Close, Yonghui Wu, Stefano Lonardi
Maintainer: Julian Taylor \<julian.taylor@adelaide.edu.au\>

## References

Wu, Y., Bhat, P., Close, T.J, Lonardi, S. (2008) Efficient and Accurate
Construction of Genetic Linkage Maps from Minimum Spanning Tree of a
Graph. Plos Genetics, **4**, Issue 10.

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## See also

[`qtl-package`](https://rdrr.io/pkg/qtl/man/a.starting.point.html)
