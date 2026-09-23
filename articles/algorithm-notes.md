# Notes on the MSTmap algorithm

> These notes refer to objects built in [Worked example
> I](https://drj001.github.io/ASMap/articles/worked-example-construction.md).
> The code below rebuilds them.

This chapter presents miscellaneous additional information associated
within the MSTmap algorithm and the arguments supplied to the R/ASMap
functions
[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
and
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md).

## MSTmap and distance calculations

After close scrutiny and experience with the MSTmap algorithm it appears
there may be circumstances where the algorithm produces inflated genetic
distances. This is especially prevalent when the linkage map being
constructed contains genotypes with many missing values. In these cases
the constructed linkage map is likely to contain runs of missing values
for these genotypes. For example, this phenomenon can be seen in the
figure below (the dark lines) for the initial constructed linkage map
`mapBC4` of the previous chapter. Seven of these genotypes were
eventually removed for having excess recombinations across the genome
but their missing value structure has also contributed substantially to
the inflated distances of the linkage groups.

``` r

plotMissing(mapBC4)
```

![Plot of the missing allele scores for the constructed map
mapBC4](algorithm-notes_files/figure-html/dist1-1.png)

Plot of the missing allele scores for the constructed map mapBC4

This problem is easier to explain by presenting various outputs of the
linkage map `mapBC4`. Firstly, the map is reconstructed but with
`return.imputed = TRUE` added to the call to ensure the imputed marker
probability matrix for the linkage groups is added to the returned cross
object.

``` r

mapBC4i <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, return.imputed = TRUE)
mapBC4i$geno[[1]]$map[1:14]
```

      mark246   mark986  mark1452  mark2266   mark391   mark426   mark670  mark1220 
    0.0000000 0.0000000 0.0000000 0.0000000 0.7641466 0.7641466 0.7641466 0.7641466 
     mark1789  mark2193  mark2269  mark2382  mark2579  mark2934 
    0.7641466 0.7641466 0.7641466 0.7641466 0.7641466 0.7641466 

``` r

mapBC4i$imputed.geno[[1]]$map[1:5]
```

      mark246   mark391   mark284   mark157  mark1990 
    0.0000000 0.7641466 1.5507414 2.3555151 3.4031908 

To exemplify the problem the markers on L.1 are used. The first 14
markers of L.1 contain the first two co-locating sets of markers and
have been calculated to be 0.76 cM apart. The imputed map for L.1 is a
reduced version of the linkage map from L.1 that contains all unique
markers and one member from each co-locating set. Markers mark246 and
mark391 have been chosen as the members for the first two co-locating
sets. It can be quickly calculated how many pairs of alleles are
observed between the two markers and how many observed recombinations
there were across these pairs.

``` r

len <- apply(mapBC4$geno[[1]]$data[,c(1,5)], 1, function(el)
             length(el[!is.na(el)]))
length(len[len > 1])
```

    [1] 281

``` r

bca <- apply(mapBC4i$geno[[1]]$data[,c(1,5)], 1, function(el){
    el <- el[!is.na(el)]
    sum(abs(diff(el)))})
bca[bca > 0]
```

    BC004 
        1 

One genuine recombination from 281 observed pairs creates a 0.35 cM
distance between the markers, falling short of 0.76 cM. The reason
behind this shortfall becomes clear when the imputed marker probability
data is printed for the first five markers of the seven genotypes with
excessive missing values

``` r

mapBC4i$imputed.geno[[1]]$data[pg$xo.lambda,1:5]
```

               mark246     mark391     mark284    mark157    mark1990
    BC009 1.000000e+00 1.000000000 1.000000000 1.00000000 1.000000000
    BC013 1.026860e-03 0.000000000 0.002317796 0.00000000 0.003386561
    BC032 1.000000e+00 0.735586234 0.672130843 0.60287457 0.639670176
    BC061 5.363308e-01 0.552323021 0.539409197 0.54183522 1.000000000
    BC069 0.000000e+00 0.008689148 0.000000000 0.02198752 0.035294621
    BC077 3.478224e-07 0.000000000 0.000000000 0.00000000 0.000000000
    BC085 9.095952e-01 0.955278155 1.000000000 0.90765213 0.881923903

Genotype BC061 has an extended set of missing values across the first
four co-located sets of markers. Using the distance formula in the
M-step of the EM algorithm given in (5), the estimates of probability
for the two missing allele scores in the first two markers for BC061
adds another 1/2 a recombination or $`~`$0.15 cM distance between the
markers. Genotype BC061 will add a similar distance between the second
and third marker set and the third and fourth marker sets due the
imputed probabilities being close to 0.5 for the run of missing values
across the four co-located sets. If there are multiple runs of missing
values for genotypes, these small distances accumulate quickly and
linkage group lengths appear inflated.

There are two solutions to this problem. The obvious first solution is
to remove the genotypes regardless off their usefulness in the genetic
map. The second less obvious solution is to use the R/ASMap function
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md) to
re-estimate the genetic distances (see section Rapid genetic distance
estimation). As the function imputes the missing values in the linkage
map with actual allele calls it circumvents the uncertain probabilities
that MSTmap uses to form its genetic distances.

``` r

mapBC4e <- quickEst(mapBC4)
chrlen(mapBC4)
```

           L.1        L.2        L.3        L.4        L.5        L.6        L.7 
    304.910957 266.240647  78.982131 252.281760  33.226962 233.485952 153.315888 
           L.8        L.9 
    106.290403   6.657526 

``` r

chrlen(mapBC4e)
```

           L.1        L.2        L.3        L.4        L.5        L.6        L.7 
    213.056491 201.340804  63.404825 187.741586  26.574889 189.370258 132.623157 
           L.8        L.9 
     77.211705   2.912879 

The difference in the lengths of the linkage groups for the two linkage
maps is dramatic with an 80+ cM reduction in L.1 and 60+ cM reductions
in L.2 and L.4.

In summary, if a linkage map is constructed with genotypes containing
many missing values then the missing value structure of the constructed
object should be checked. If runs of missing values are detected then
the lengths of linkage groups are most likely inflated and estimation of
genetic distances should be checked using
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md).

## Use of the argument `mvest.bc`

The inflated distances caused by excessive missing values across the
marker set for a subset of genotypes can be further exacerbated by an
injudicious use of the argument `mvest.bc` in the R/ASMap construction
functions
[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
and
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md).

To understand how this is possible, there is a requirement to understand
how the MSTmap algorithm determines co-segregating or co-locating
markers. In the initial stages of the algorithm two markers are deemed
co-located if the pairwise distance is zero with this calculation
occurring only from pairs of alleles that are observed in both markers.
For each co-locating set a marker is chosen as the representative marker
and the remaining markers in the set are placed aside with a knowledge
of their direct link to the representative marker. The linkage map is
then constructed with only the representative markers. This scenario is
equivalent to setting `mvest.bc = FALSE`.

If `mvest.bc = TRUE` then missing allele calls contained in each marker
are estimated to have probability 0.5 of being an A allele *before*
clustering of the markers has occurred. This has an obvious advantage of
involving the complete set of genotypes when calculating pairwise
information between markers. However, for any two markers where there is
some proportion of unobserved allele pairs, the pairwise distance
between the markers becomes non-zero and they are deemed to be not
co-located. This results in a larger set of markers being involved in
the linkage map construction with the original missing allele calls
being continually imputed as part of the marker ordering algorithm (see
section Marker Ordering). For situations where there is an increased
number of missing values across a set of markers for a subset of
genotypes, the cumulative distances between adjacent markers quickly
increases.

To illustrate this problem, the barley backcross marker set `mapBC3` is
constructed using the function
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
with the argument `mvest.bc = TRUE`.

``` r

mapBC4a <- mstmap(mapBC3, bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, mvest.bc = TRUE)
nmar(mapBC4)
```

    L.1 L.2 L.3 L.4 L.5 L.6 L.7 L.8 L.9 
    574 473  75 551  32 188 160 110  10 

The clustering algorithm produces an identical number of linkage groups
with the same markers within each linkage group. However, the number of
uniquely located markers within each linkage group increases
dramatically between `mapBC4` and `mapBC4a`. This increase in the number
of uniquely located markers also increases the length of each linkage
group.

``` r

sapply(mapBC4a$geno, function(el) length(unique(round(el$map, 4)))) - sapply(mapBC4$geno, function(el) length(unique(round(el$map, 4))))
```

    L.1 L.2 L.3 L.4 L.5 L.6 L.7 L.8 L.9 
    384 208  44 305  14  97  62  49   4 

``` r

chrlen(mapBC4a)
```

           L.1        L.2        L.3        L.4        L.5        L.6        L.7 
    498.063994 340.716640  96.513164 356.269346  41.016980 275.675999 164.254019 
           L.8        L.9 
    123.280455   6.952472 

Again, the solution to this problem is to use the efficient distance
estimation function
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md) to
provide an estimate of the genetic distances of the markers in each of
the linkage groups.

``` r

mapBC4b <- quickEst(mapBC4a)
chrlen(mapBC4b)
```

          L.1       L.2       L.3       L.4       L.5       L.6       L.7       L.8 
    207.88659 199.40969  62.42760 184.52014  25.92514 187.42973 132.71937  75.91499 
          L.9 
      2.58925 

## Use of the argument `detectBadData`

The last two sections discussed constructed linkage map scenarios where
the genetic distance calculations became inflated due to an excessive
number of missing values being present for a subset of genotypes in the
marker set. If this situation is not present and the linkage group
distances are inflated then it is most likely caused by an increased
recombination rate in some or all of the genotypes. This is easily
checked using the appropriate call to
[`statGen()`](https://drj001.github.io/ASMap/reference/statGen.md) or
[`profileGen()`](https://drj001.github.io/ASMap/reference/profileGen.md)
after an initial linkage map is constructed. If there is an increased
recombination rate in a small subset of genotypes they can be removed
before further linkage map construction. If there appears to be a
general inflation in the recombination rate across a large set or all of
the genotypes then there may be genotyping errors that have occurred in
the physical process used to obtain the allele calls for the markers.
Consequently, reduction of the linkage group lengths becomes more
problematic.

To circumvent this issue, the construction functions
[`mstmap.data.frame()`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
and
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
in the R/ASMap package contain an argument called `detectBadData` that,
if set to `TRUE`, instructs MSTmap to detect suspiciously called alleles
(see equation 6 of section Marker Ordering for more details). If
suspicious allele calls are found they are set to missing and imputed
using the EM algorithm detailed in section Marker Ordering. This creates
a reduction in the number of recombinations between adjacent markers.
The consequence of this is a possibly dramatic reduction in the
distances between adjacent markers and overall linkage group lengths.

To exemplify the use of `detectBadData`, simulated genotyping errors are
added to the final barley backcross linkage map `mapBC` by randomly
switching singular allele calls in the data for each linkage group.

``` r

mapBCd <- mapBC
mapBCd$geno <- lapply(mapBCd$geno, function(el){
    ns <- sample(1:ncol(el$data), ncol(el$data)/2, replace = TRUE)
    ns <- cbind(sample(1:nrow(el$data), ncol(el$data)/2, replace = TRUE), ns)
    el$data[ns] <- abs(1 - el$data[ns])
    el$data[el$data == 0] <- 2
    el})
mapBCd <- quickEst(mapBCd)
chrlen(mapBCd)
```

         L.1      L.2      L.3      L.4      L.5      L.6      L.7 
    262.7372 265.0427 169.6948 252.3984 223.7460 168.1366 153.2918 

The function
[`quickEst()`](https://drj001.github.io/ASMap/reference/quickEst.md)
will not ignore the genotyping errors that have been added to the data
and linkage groups will appear inflated.

``` r

mapBCda <- mstmap(mapBCd, bychr = TRUE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-12, detectBadData = TRUE)
chrlen(mapBCda)
```

         L.1      L.2      L.3      L.4      L.5      L.6      L.7 
    228.9131 224.6887 141.4048 205.6216 196.4040 142.0062 137.3560 

A judicious use of the `detectBadData = TRUE` in the
[`mstmap.cross()`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)
call instructs the algorithm to detect and ignore these errors during
the marker ordering algorithm of MSTmap. This deflates the linkage group
distances to similar lengths as the linkage groups lengths of the final
linkage map `mapBC`.
