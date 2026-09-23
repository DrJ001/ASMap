# How the MSTmap algorithm works

## Background

Genetic linkage maps are widely used in the biological research
community to explore the underlying DNA of populations. They generally
consist of a set of polymorphic genetic markers spanning the entire
genome of a population generated from a specific cross of parental
lines. This exploration may involve the dissection of the linkage map
itself to understand the genetic landscape of the population or, more
commonly, it is used to conduct gene-trait associations such as
quantitative trait loci (QTL) analysis or genomic selection (GS). For
QTL analysis, the interpretation of significant genomic locations is
enhanced if the linkage map contains markers that have been assigned and
optimally ordered within chromosomal groups. This can be achieved
algorithmically by using linkage map construction techniques that
utilise assumptions of Mendelian genetics.

In the R statistical computing environment ([R Development Core Team
2014](#ref-rsoft14)) there is only a handful of packages that can
perform linkage map construction. A popular package is the linkage map
construction and QTL analysis package R/qtl ([Broman and Wu
2014](#ref-br14)). Since its inception in late 2001 it has grown
considerably and incorporates functionality for a wide set of
populations. These include a simple Backcross (BC), Doubled Haploid
(DH), intercrossed F2 (F2) and 4-way crosses. More recently, the authors
have added functionality for additional populations generated from
advanced Recombinant Inbred Lines (RIL) as well as populations generated
from multiple backcrossing and selfing processes. The package offers an
almost complete set of tools to construct, explore and manipulate
genetic linkage maps for each population. To complement the package the
authors have also publihed a comptrehensive book ([Broman and Sen
2009](#ref-brosen09)). A more recent available addition to the
contributed package list in R is R/onemap ([Margarido and Mollinari
2014](#ref-one14)). This package provides a suite of tools for genetic
linkage map construction for BC, DH and RIL inbred populations and also
has functionality for outcrossing populations.

Both packages contain functions that perform moderately well at
clustering markers into homogeneous linkage groups. Unfortunately, the
functions in R/qtl and R/onemap used to perform optimal marker ordering
within linkage groups are tediously slow. In both packages there are
functions that compute initial orders using non-exhaustive methods such
as SERIATION ([Buetow and Chakravarti 1987](#ref-ser87)) and RECORD
([Van Os et al. 2005](#ref-rec05)). As these methods may not return an
optimal order, a second stage antiquated brute force approach is used
that checks all order combinations within a marker window. This is known
throughout the linkage map construction community as “rippling”. For
linkage groups with large numbers of markers and a moderate sized window
(eight markers), rippling can be very computationally cumbersome and
also may need to be performed several times before an optimal order is
reached.

In an attempt to circumvent these computational issues the R/ASMap
package was developed ([Taylor and Butler 2014](#ref-tb14)). The package
contains linkage map construction functions that utilise the MSTmap
algorithm derived in ([Wu 2008](#ref-mst08)) and computationally
implemented as C++ source code freely available at
[http://alumni.cs.ucr.edu/~yonghui/mstmap.html](http://alumni.cs.ucr.edu/~yonghui/mstmap.md).
The algorithm uses the minimum spanning tree of a graph ([Cheriton and
Tarjan 1976](#ref-ct76)) to cluster markers into linkage groups as well
as find the optimal marker order within each linkage group in a very
computationally efficient manner. In contrast to R/qtl and R/onemap,
genetic linkage maps are constructed using a one-stage approach. The
algorithm is restricted to linkage map construction with Backcross (BC),
Doubled Haploid (DH) and Recombinant Inbred (RIL) populations. For RIL
populations, the level of self pollination can also be given and
consequently the algorithm can handle selfed F2, F3, $`\ldots`$, Fn
populations where $`n`$ is the level of selfing. Advanced RIL
populations are also allowed as they are treated like a bi-parental for
the purpose of linkage map construction.

The R/ASMap package uses an R/qtl format for the structure of its
genetic objects. Once the class of the object is appropriately set, both
R/ASMap and R/qtl functions can be used synergistically to construct,
explore and manipulate the object. To complement the efficient MSTmap
linkage map construction functions the R/ASMap package also contains a
function that “pulls” markers of different types from the linkage map
and places them aside. A complementary “push” function also exists to
push the markers back to linkage groups ready for construction or
reconstruction. There are also several numerical and graphical
diagnostic tools to efficiently check the quality of the constructed
map. This includes the ability to simultaneously graphically display
multiple panel profiles of linkage map statistics for the set of
genotypes used. Additionally, profile marker/interval statistics can be
simultaneously displayed for the entire genome or subsetted to
pre-defined linkage groups. In tandem, these fast graphical diagnostic
tools and efficient linkage map construction assist in providing a rapid
turnaround time in linkage map construction and diagnosis.

The vignette is set out in chapters. The final section of the
introduction provides a brief technical explanation of the components of
the MSTmap algorithm. In the second chapter the functions of R/ASMap are
discussed in detail and examples are provided, where appropriate, using
a data set integrated into the package. Chapter 3 contains a completely
worked example of a barley Backcross population including
pre-construction diagnostics, linkage map construction, and
post-construction diagnostics using the available R/ASMap functions. The
chapter also discusses how R/ASMap can be used for post construction
linkage map improvement through techniques such as fine mapping or
combining linkage maps. The last chapter presents some additional useful
information on aspects of the MSTmap algorithm learnt through detailed
exploration of the linkage map construction functions in this package.

## MSTmap

It is important to understand some of the technical features of the
MSTmap algorithm as they are contained in the two linkage map
construction functions that are available with the R/ASMap package.

Following the notation of ([Wu 2008](#ref-mst08)) consider a Doubled
Haploid population of $`n`$ individuals genotyped across a set of $`t`$
markers where each $`(i,j)`$th entry of the $`n \times t`$ matrix
$`\boldsymbol{M}`$ is either an $`A`$ or a $`B`$ representing the two
parental homozygotes in the population. Let $`\boldsymbol{P}_{jk}`$ be
the probability of a recombination event between the markers
$`(\boldsymbol{m}_j,
\boldsymbol{m}_k)`$ where $`0 \leq \boldsymbol{P}_{jk} < 0.5`$. MSTmap
uses two possible weight objective functions based on recombination
probabilities between the markers
``` math
\begin{aligned}
  w_p(j,k) &= \boldsymbol{P}_{jk} \qquad (1)\\
  w_{ml}(j, k) &= -\bigl(\boldsymbol{P}_{jk}\log\,\boldsymbol{P}_{jk} + (1 -
  \boldsymbol{P}_{jk})\log(\,1 - \boldsymbol{P}_{jk})\bigr) \qquad (2)
\end{aligned}
```
In general $`\boldsymbol{P}_{jk}`$ is not known and so it is replaced by
an estimate, $`d_{jk}/n`$ where $`d_{jk}`$ corresponds to the hamming
distance between $`\boldsymbol{m}_j`$ and $`\boldsymbol{m}_k`$ (the
number of non-mathcing alleles between the two markers). This estimate,
$`d_{jk}/n`$, is also the maximum likelihood estimate for
$`\boldsymbol{P}_{jk}`$ for the two weight functions defined above.

### Clustering

If markers $`\boldsymbol{m}_j`$ and $`\boldsymbol{m}_k`$ belong to two
different linkage groups then $`\boldsymbol{P}_{jk} = 0.5`$ and the
hamming distance between them has the property $`E(d_{jk}) = n/2`$. A
simple thresholding mechanism to determine whether markers belong to the
same linkage group can be calculated using Hoeffdings inequality, \$\$
``` math
\begin{aligned}
 P(d_{jk} < \delta) \leq \mbox{exp}(-2(n/2 - \delta)^2/n)
 \qquad (3)
 
\end{aligned}
```

\$\$ for $`\delta < 0.5`$. For a given $`P(d_{jk} < \delta) = \epsilon`$
and $`n`$, the equation $`-2(n/2 - \delta)^2/n = \mbox{log}\,\epsilon`$
is solved to determine an appropriate hamming distance threshold,
$`\hat{\delta}`$. ([Wu 2008](#ref-mst08)) indicate that the choice of
$`\epsilon`$ is not crucial when attempting to form linkage groups.
However, the equation that requires solving is highly dependent on the
number of individuals in the population. For example, for a DH
populationm, the figure below shows the profiles of the
$`-\log 10 \epsilon`$ against the number of individuals in the
population for four threshold minimum cM distances $`(25, 30, 35,
40)`$. MSTmap uses a default of $`\epsilon = 0.00001`$ which would work
universally well for population sizes of
$`n \sim 150\, \mbox{to}\, 200`$. For larger numbers of individuals, for
example 350, the plot indicates an $`\epsilon = 1.0e^{-12}`$ to
$`1.0e^{-15}`$ would use a conservative minimum threshold of 30-35 cM
before linking markers between clusters. If the default $`\epsilon =
0.00001`$ is given in this instance this threshold is dropped to
$`\sim 45`$ cM and consequently distinct clusters of markers will appear
linked. For this reason, the figure below should always initially be
checked before linkage map construction to ensure an appropriate p-value
is given to the MSTmap algorithm.

To cluster the markers MSTmap uses an edge-weighted undirected complete
graph, for $`\boldsymbol{M}`$ where the individual markers are vertices
and the edges between any two markers $`\boldsymbol{m}_j`$ and
$`\boldsymbol{m}_k`$ is the pairwise hamming distance $`d_{jk}`$. Edges
with weights greater than $`\hat{\delta}`$ are then removed. The
remaining connected components allow the marker set $`\boldsymbol{M}`$
to be partitioned into $`r`$ linkage groups, $`\boldsymbol{M} =
[\boldsymbol{M}_1,\ldots,\boldsymbol{M}_r]`$.

``` r

pValue(dist = c(25, 30, 35, 40), pop.size = 100:500, map.function = "kosambi")
```

![Negative log10 epsilon versus the number of genotypes in the
population for four threshold cM
distances.](mstmap-algorithm_files/figure-html/pvalfig-1.png)

Negative log10 epsilon versus the number of genotypes in the population
for four threshold cM distances.

### Marker Ordering

For simplicity, consider the $`n \times t`$ matrix of markers
$`\boldsymbol{M}`$ belongs to the same linkage group. Preceding marker
ordering, the markers are “binned” into groups where, within each group,
the pairwise distance between any two markers is zero. The markers
within each group have no recombinations between them and are said to be
co-locating at the same genomic location for the $`n`$ genotypes used to
construct the linkage map. A representative marker is then chosen from
each of the bins and used to form the reduced $`n \times t^*`$ marker
set $`\boldsymbol{M}^*`$.

For the reduced matrix $`\boldsymbol{M}^*`$, consider the complete set
of entries $`(j, k)
\in (1, \ldots, t^*)`$ for either weight function (1) or (2). These
complete set of entries can be viewed as the upper triangle of a
symmetric weight matrix $`\boldsymbol{W}`$. MSTmap views all these
entries as being connected edges in an undirected graph where the
individual markers are vertices. A marker order for the set
$`\boldsymbol{M}^*`$, also known as a travelling salesman path (TSP),
can be determined by visiting each marker once and summing the weights
from the connected edges. To find a minimum weight (TSP$`_{min}`$),
MSTmap uses a minimum spanning tree (MST) algorithm ([Cheriton and
Tarjan 1976](#ref-ct76)), such as Prims algorithm ([Prim
1957](#ref-prim57)). If the TSP$`_{min}`$ is unique then the MST is the
correct order for the markers. For cases where the data contains
genotyping errors or lower numbers of individuals the MST may not be a
complete path and contain markers or small sets of markers as individual
nodes connected to the path. In these cases, MSTmap uses the longest
path in the MST as the backbone and employs several efficient local
optimization techniques such as K-opt, node-relocation and
block-optimize (see [Wu 2008](#ref-mst08)) to improve the current
minimum TSP. By integrating these local optimization techniques into the
algorithm, MSTmap provides users with a true one stage marker ordering
algorithm.

One excellent feature of the MSTmap algorithm is the utilisation of an
EM type algorithm for the imputation of missing allele scores that is
tightly integrated with the ordering algorithm for the markers. To
achieve this the marker matrix $`\boldsymbol{M}^*`$ is converted to a
matrix, $`\boldsymbol{A}`$, where the entries represent the
probabilistic certainty of the allele being A. For the $`j`$th marker
and $`i`$th individual then
``` math
\begin{aligned}
  \boldsymbol{A}(i,j) = \left\{\begin{array}{ll}
          1 & \quad \mbox{if $\boldsymbol{M}^*(i,j)$ is the A allele}\\
          0 & \quad \mbox{if $\boldsymbol{M}^*(i,j)$ is the B allele}\\
          \frac{(1 - \hat{\boldsymbol{P}}_{j - 1,j})(1 - \hat{\boldsymbol{P}}_{j,j+1})}
          {(1 - \hat{\boldsymbol{P}}_{j,j+1})(1 - \hat{\boldsymbol{P}}_{j,j+1})
           + \hat{\boldsymbol{P}}_{j-1,j} \hat{\boldsymbol{P}}_{j,j+1}} & \quad \mbox{if
           $\boldsymbol{M}^*(i,j)$ is missing} \end{array}\right.\qquad (4)
\end{aligned}
```
where $`\hat{\boldsymbol{P}}_{j,j-1}`$ and
$`\hat{\boldsymbol{P}}_{j,j+1}`$ are estimated recombination fractions
between the $`(j-1)`$th and $`j`$th marker and $`j`$th and $`(j+1)`$th
marker respectively. The equation on the right hand side is the
posterior probability of the missing value in marker $`j`$ being the A
allele for genotype $`i`$ given the current estimate. The ordering
algorithm begins by initially calculating pairwise normalized distances
between all markers in $`\boldsymbol{M}^*`$ and deriving an initial
weight matrix, $`\boldsymbol{W}`$. An undirected graph is formed using
the markers as vertices and the upper triangular entries of
$`\boldsymbol{W}`$ as connected edges. An MST of the undirected graph is
then found to establish an initial order for the markers of the linkage
group. For the current order at the $`(j - 1, j,
j + 1)`$th markers the E-step of algorithm requires updating the missing
observation at marker $`j`$ by updating the estimates
$`\hat{\boldsymbol{P}}_{j-1,j} = \hat{d}_{j-1,j}/n`$ and
$`\hat{\boldsymbol{P}}_{j,j+1} =
    \hat{d}_{j,j+1}/n`$ in (4). The M-step then re-estimates the
pairwise distances between all markers in $`\boldsymbol{M}^*`$ where,
for the $`j`$th and $`k`$th marker, is
``` math
\begin{aligned}
\hat{d}_{jk} = \sum_{i = 1}^{t^*}\boldsymbol{A}(i,j)(1 - \boldsymbol{A}(i, k)) + \boldsymbol{A}(i,k)(1 - \boldsymbol{A}(i, j))
\qquad (5)
\end{aligned}
```
and the weight matrix $`\boldsymbol{W}`$ is recalculated. An undirected
graph is formed with the markers as vertices and the upper triangular
entries of $`\boldsymbol{W}`$ as connected edges. A new order of the
markers is derived by obtaining an MST of the undirected graph and the
algorithm is repeated to convergence. Although this requires several
iterations to converge, the computational time for the ordering
algorithm remains expedient. However, an increase in the number of
missing values will increase computation time.

If required, the MSTmap algorithm also detects and removes genotyping
errors as well as integrates this process into the ordering algorithm.
The technique involves using a weighted average of nearby markers to
determine the expected state of the allele. For individual $`i`$ and
marker $`j`$ the expected value of the allele is calculated using
``` math
\begin{aligned}
\mbox{E}[\boldsymbol{A}(i, j)] = \sum_{j \neq
    k}d_{j,k}^{-2}\boldsymbol{A}(i,k)\bigg/\sum_{j \neq k}d_{j,k}^{-2}
\qquad (6)
\end{aligned}
```
In this equation the weights are the inverse square of the distance from
marker $`j`$ to its nearby markers. MSTmap only uses a small set of
nearby markers during each iteration and the observed allele is
considered suspicious if
$`|\mbox{E}[\boldsymbol{A}(i, j)] - \boldsymbol{A}(i, j)| >
0.75`$. If an observation is detected as suspicious it is treated as
missing and imputed using the EM algorithm discussed previously. The
removal of the suspicious allele has the effect of reducing the number
of recombinations between the marker containing the suspicious
observation and the neighbouring markers. This has an influential effect
on the genetic distance between markers and the overall length of the
linkage group.

The complete algorithm used to initially cluster the markers into
linkage groups and optimally order markers within each linkage group,
including imputing missing alleles and error detection, is known as the
MSTmap algorithm.

Broman, K. W, and S Sen. 2009. *A Guide to QTL Mapping with /*.
Springer-Verlag.

Broman, K. W, and H Wu. 2014. *: Tools for Analayzing QTL Experiments*.
<http://www.CRAN.R-project.org/src/contrib/Archive/qtl/>.

Buetow, K. H, and A Chakravarti. 1987. “Multipoint gene mapping using
seriation. 1 Gneeral methods.” *American Journal of Human Genetics* 41:
180–88.

Cheriton, D., and R. E. Tarjan. 1976. “Finding Minimum Spanning Trees.”
*SIAM Journal on Computing* 5 (4): 724–42.
<https://doi.org/10.1137/0205051>.

Margarido, G, and M Mollinari. 2014. *: Software for constructing
genetic maps in experimental crosses: full-sib, RILs, F2 and
backcrosses*. <http://www.CRAN.R-project.org/src/contrib/Archive/qtl/>.

Prim, R. C. 1957. “Shortest Connection Networks and Some
Generalizations.” *Bell System Techical. Journal* 36: 1389–401.

R Development Core Team. 2014. *R: A Language and Environment for
Statistical Computing*. R Foundation for Statistical Computing.
<http://www.R-project.org>.

Taylor, J. D, and D Butler. 2014. *: An (A)ccurate and (S)peedy linkage
map construction package for inbred populatons that uses the extremely
efficient MSTmap algorithm.*
<http://www.CRAN.R-project.org/src/contrib/Archive/ASMap/>.

Van Os, Hans, Piet Stam, Richard Visser, and Herman Van Eck. 2005.
“RECORD: A Novel Method for Ordering Loci on a Genetic Linkage Map.”
*TAG Theoretical and Applied Genetics* 112: 30–40.
<http://dx.doi.org/10.1007/s00122-005-0097-x>.

Wu, Prasanna R. AND Close, Yonghui AND Bhat. 2008. “Efficient and
Accurate Construction of Genetic Linkage Maps from the Minimum Spanning
Tree of a Graph.” *PLoS Genetics* 4 (10).
<https://doi.org/10.1371/journal.pgen.1000212>.
