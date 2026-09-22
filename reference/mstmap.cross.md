# Extremely fast linkage map construction for qtl objects using MSTmap.

Extremely fast linkage map construction for qtl objects using the source
code for MSTmap (see Wu et al., 2008). The construction includes linkage
group clustering, marker ordering and genetic distance calculations.

## Usage

``` r
# S3 method for class 'cross'
mstmap(object, chr, id = "Genotype", bychr = TRUE,
       suffix = "numeric", anchor = FALSE, dist.fun = "kosambi",
       objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
       noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE,
       detectBadData = FALSE, return.imputed = FALSE,
       trace = FALSE, ...)
```

## Arguments

- object:

  A `"cross"` object generated from the qtl package. Specifically the
  object needs to inherit from one of the following classes `"bc"`,
  `"dh"`, `"riself"`, `"bcsft"` (see Details).

- chr:

  A character string of linkage group names that require re-construction
  and/or optimal ordering of the markers they contain. (see Details).

- id:

  The name of the column in `object$pheno` that uniquely identifies the
  genotype names. Default is `"Genotype"`.

- bychr:

  Logical value. For a given set of linkage groups defined by `chr`, if
  `TRUE` then split linkage groups (only if required, see `p.value`) and
  order markers within linkage groups. If `FALSE` then combine linkage
  groups and reconstruct. Default is `TRUE`.

- suffix:

  Character string either `"numeric"` or `"alpha"` determining whether
  numeric or alphabetic ascending values are post-fixed to linkage group
  names when splitting linkage groups.

- anchor:

  Logical value. The MSTmap algorithm does not respect the inputted
  marker order of the linkage map required for construction. For a given
  set of linkage groups defined by `chr`, if `TRUE` the order of the
  inputted markers is respected regardless of the choices of `chr` and
  `bychr`. Default is `FALSE`.

- dist.fun:

  Character string defining the distance function used for calculation
  of genetic distances. Options are "kosambi" and "haldane". Default is
  "kosambi".

- objective.fun:

  Character string defining the objective function to be used when
  constructing the map. Options are `"COUNT"` for minimising the sum of
  recombination events between markers and `"ML"` for maximising the
  likelihood objective function. Default is `"COUNT"`.

- p.value:

  Numerical value to specify the threshold to use when clustering
  markers. Defaults to `1e-06`. If a value greater than one is given
  this feature is turned off inputted marker data are assumed to belong
  to the same linkage group (see Details).

- noMap.dist:

  Numerical value to specify the smallest genetic distance a set of
  isolated markers can appear distinct from other linked markers.
  Isolated markers will appear in their own linkage groups and will be
  of size specified by `noMap.size`.

- noMap.size:

  Numerical value to specify the maximum size of isolated marker linkage
  groups that have been identified using `noMap.dist`. This feature can
  be turned off by setting it to 0. Default is 0.

- miss.thresh:

  Numerical value to specify the threshold proportion of missing marker
  scores allowable in each of the markers. Markers above this threshold
  will not be included in the linkage map. Default is 1.

- mvest.bc:

  Logical value. If `TRUE` missing markers will be imputed before
  clustering the markers into linkage groups. This is restricted to
  `"bc","dh","riself"` populations only (see Details). Default is
  `FALSE`.

- detectBadData:

  Logical value. If `TRUE` possible genotyping errors are detected, set
  to missing and then imputed as part of the marker ordering algorithm.
  Genotyping errors will also be printed in the file specified by
  `trace`. This is restricted to `"bc","dh","riself"` populations only.
  (see Details). Default is `FALSE`.

- return.imputed:

  Logical value. If `TRUE` then the imputed marker probability matrix is
  returned for the linkage groups that are constructed (see Details).
  Default is `FALSE`.

- trace:

  An automatic tracing facility. If `trace = FALSE` then minimal
  `MSTmap` output is piped to the screen during the algorithm. If
  `trace = TRUE`, then detailed output from MSTmap is piped to
  `"MSToutput.txt"`. This file is equivalent to the output that would be
  obtained from running the MSTmap executable from the command line.

- ...:

  Currently ignored.

## Details

The qtl cross object needs to inherit one of the allowable classes
`"bc","dh","riself", "bcsft"`. This provides a safeguard against
attempts to construct a map for more complex populations that can exist
in qtl. Users should be aware when doubled haploid populations are read
in using [`read.cross()`](https://rdrr.io/pkg/qtl/man/read.cross.html)
from the qtl package they inherit the class `"bc"`. Users can apply the
class `"dh"` by simply changing the class of the object. For the purpose
of linkage map construction the classes `"bc"` and `"dh"` will provide
equivalent results.

MSTmap supports `"RILn"` populations, where n is the number of
generations of selfing. Markers in these populations are required to be
fully informative i.e. contain 3 distinct allele types such as AA, BB
for parental homozygotes and AB for phase unknown heterozygotes. If
`read.cross` is used to import the `"RILn"` population the resultant
object will initially be given a class `"f2"`. The level of selfing
would then have to be encoded into the object by applying one of the two
conversion functions available in the qtl package. For a population that
has been generated by selfing n times the conversion function
`convert2bcsft` can be used by setting the arguments `F.gen = n` and
`BC.gen = 0`. Populations that are genuine advanced RILs can be
converted using the `convert2riself` function.

This method function is designed to be an "all-in-one" function that
will allow you to construct linkage maps extremely fast in multiple
different ways from the supplied cross `object`. Initially, the map can
be kept complete or a subset of selected linkage groups can be chosen
using the `chr` argument. Setting `bychr = FALSE` will bulk the marker
information for the selected linkage groups and, if necessary, form new
linkage groups and optimise the marker order within each. Setting
`bychr = TRUE` will ensure that markers are optimally ordered within
each linkage group. This will also break linkage groups depending on the
p-value given in the call (see below for details of the use of
`p.value`). If the linkage map was initially subsetted, the linkage
groups not involved in the subset are returned to ensure the map is
complete.

The algorithm allows an adjustment of the `p.value` threshold for
clustering of markers to distinct linkage groups (see Wu et al., 2008)
and is highly dependent on the number of individuals in the population.
As the number of individuals increases the `p.value` threshold should be
decreased accordingly. This may require some trial and error to achieve
desired results. When `bychr = TRUE`, established linkage groups may
also split depending on the `p.value` given. To prevent this the
`p.value` threshold may be increased to a desired value or the splitting
may be prevented altogether by supplying a value greater than one to
this argument.

If `mvest.bc = TRUE` and the population type is `"bc","dh","riself"`
then missing values are imputed before markers are clustered into
linkage groups. This is only a simple imputation that places a 0.5
probability of the missing observation being one allele or the other and
is used to assist the clustering algorithm when there is known to be
high numbers of missing observations between pairs of markers.

It should be highlighted that for population types `"bc","dh","riself"`,
imputation of missing values occurs regardless of the value of
`mvest.bc`. This is achieved using an EM algorithm that is tightly
coupled with marker ordering (see Wu et al., 2008). Initially a marker
order is obtained omitting missing marker scores and then imputation is
performed based on the underlying recombinant probabilities of the
flanking markers with the markers containing the missing value. The
recombinant probabilities are then recomputed and an update of the
pairwise distances are calculated. The ordering algorithm is then run
again and the complete process is repeated until convergence. Note, the
imputed probability matrix for the linkage map being constructed is
returned if `return.imputed = TRUE`.

For populations `"bc","dh","riself"`, if `detectBadData = TRUE` the
marker ordering algorithm also includes the detection of genotyping
errors. For any individual genotype, the detection method is based on a
weighted Euclidean metric (see Wu et al., 2008) that is a function of
the recombination probabilities of all the markers with the marker
containing the suspicious observation. Any genotyping errors detected
are set to missing and the missing values are then imputed as part of
the marker ordering algorithm. Note, the detection of these errors and
their amendment can be returned in the imputed probability matrix if
`return.imputed = TRUE`.

If `return.imputed = TRUE` and the object has class `"bc","dh","riself"`
then the marker probability matrix is returned for the linkage groups
that have been constructed using the algorithm. Each linkage group is
named identically to the linkage groups of the map and contains an
ordered `"map"` element and a `"data"` element consisting of marker
probabilities of the A allele being present (i.e. P(A) = 1, P(B) = 0).
Both elements contain a possibly reduced version of the marker set that
includes all non-colocating markers as well as the first marker of any
set of co-locating markers.

## Value

The function returns a cross object with an identical class structure to
the cross `object` inputted. The object is a list with usual components
`"pheno"` and `"geno"`. If markers were omitted for any reason during
the construction, the object will have an `"omit"` component with all
omitted markers in a collated matrix. If `return.imputed = TRUE` then
the object will also contain an `"imputed.geno"` element.

## References

Wu, Y., Bhat, P., Close, T.J, Lonardi, S. (2008) Efficient and Accurate
Construction of Genetic Linkage Maps from Minimum Spanning Tree of a
Graph. Plos Genetics, **4**, Issue 10.

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor, Dave Butler, Timothy Close, Yonghui Wu, Stefano Lonardi

## See also

[`mstmap.data.frame`](https://drj001.github.io/ASMap/reference/mstmap.data.frame.md)
and
[`breakCross`](https://drj001.github.io/ASMap/reference/breakCross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## bulking linkage groups and reconstructing entire linkage map

test1 <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi",
                trace = FALSE)
#> Number of linkage groups: 24
#> The size of the linkage groups are: 41   5   33  10  40  35  8   13  37  30  6   9   21  32  6   54  56  6   27  41  15  33  37  4   
#> The number of bins in each linkage group: 41 4   33  9   40  34  7   13  36  29  5   8   21  32  6   53  55  5   27  41  14  32  36  3   
pull.map(test1)
#> $L.1
#>     1A.m.41     1A.m.40     1A.m.39     1A.m.38     1A.m.37     1A.m.36 
#>   0.0000000   0.4587285   0.9174569   3.2126457  10.1373152  11.5138097 
#>     1A.m.35     1A.m.34     1A.m.33     1A.m.32     1A.m.31     1A.m.30 
#>  16.1139008  23.9762340  24.4349625  24.8936909  25.9826349  27.1876743 
#>     1A.m.29     1A.m.28     1A.m.27     1A.m.26     1A.m.25     1A.m.24 
#>  29.0233610  32.2387955  32.6975239  34.5332107  63.5518124  64.9283068 
#>     1A.m.23     1A.m.22     1A.m.21     1A.m.20     1A.m.19     1A.m.18 
#>  68.1437413  68.6024697  69.0611982  69.9787324  71.3552268  71.8139553 
#>     1A.m.17     1A.m.16     1A.m.15     1A.m.13     1A.m.14     1A.m.12 
#>  72.7314894  73.1902179  74.5667123  75.5029664  75.9804099  76.4391384 
#>     1A.m.11     1A.m.10      1A.m.9      1A.m.8      1A.m.7      1A.m.6 
#>  77.8156328  78.7331670  81.4882455  81.9469740  82.4057024  82.8644309 
#>      1A.m.5      1A.m.4      1A.m.3      1A.m.2      1A.m.1 
#>  83.3231594  84.6996538  99.8235155 102.5785940 103.9550885 
#> 
#> $L.10
#>   3B.m.30   3B.m.29   3B.m.28   3B.m.27   3B.m.26   3B.m.25   3B.m.24   3B.m.23 
#>   0.00000  10.23223  22.39308  30.25542  34.39367  34.85278  38.99064  53.45212 
#>   3B.m.22   3B.m.21   3B.m.20   3B.m.19   3B.m.17   3B.m.18   3B.m.16   3B.m.15 
#>  56.35658  61.41969  66.48280  67.40033  68.31786  68.77659  69.23532  69.69405 
#>   3B.m.14   3B.m.13   3B.m.11   3B.m.12   3B.m.10    3B.m.9    3B.m.8    3B.m.7 
#>  73.83191  74.29064  74.74937  74.74937  75.20810  97.15983  98.07737 104.62584 
#>    3B.m.6    3B.m.5    3B.m.4    3B.m.3    3B.m.2    3B.m.1 
#> 129.29707 131.59226 132.05099 132.50972 133.42725 133.88598 
#> 
#> $L.11
#>    3D.m.6    3D.m.5    3D.m.4    3D.m.3    3D.m.1    3D.m.2 
#> 0.0000000 0.4587314 5.0589782 8.2768742 8.7379031 8.7379031 
#> 
#> $L.12
#>   4A.m.9   4A.m.7   4A.m.8   4A.m.6   4A.m.5   4A.m.3   4A.m.4   4A.m.2 
#> 0.000000 1.835687 1.835687 2.294415 2.753144 3.670678 4.129406 6.424595 
#>   4A.m.1 
#> 7.801089 
#> 
#> $L.13
#>  4A.m.30  4A.m.29  4A.m.28  4A.m.27  4A.m.26  4A.m.25  4A.m.24  4A.m.23 
#>  0.00000  6.92467 18.60024 26.93343 43.06735 44.44384 44.90257 45.82010 
#>  4A.m.22  4A.m.21  4A.m.20  4A.m.19  4A.m.18  4A.m.17  4A.m.16  4A.m.15 
#> 46.73764 47.19637 48.11390 51.79024 53.16673 53.62546 83.92398 84.38271 
#>  4A.m.14  4A.m.13  4A.m.12  4A.m.11  4A.m.10 
#> 86.67790 94.54023 94.99896 95.91650 98.21169 
#> 
#> $L.14
#>    4B.m.32    4B.m.31    4B.m.30    4B.m.29    4B.m.28    4B.m.27    4B.m.26 
#>   0.000000   8.333189  21.963866  22.881400  24.257894  29.784884  37.647217 
#>    4B.m.25    4B.m.24    4B.m.23    4B.m.22    4B.m.21    4B.m.20    4B.m.19 
#>  41.785078  42.243806  43.161340  50.086010  51.003544  52.380039  52.838767 
#>    4B.m.18    4B.m.17    4B.m.16    4B.m.15    4B.m.14    4B.m.13    4B.m.12 
#>  53.756301  56.051490  56.510218  58.345905  58.804634  59.263362  59.722091 
#>    4B.m.11    4B.m.10     4B.m.9     4B.m.8     4B.m.7     4B.m.6     4B.m.5 
#>  60.180819  60.639547  61.557082  62.015810  64.310999  78.436651  81.191729 
#>     4B.m.4     4B.m.3     4B.m.2     4B.m.1 
#>  83.946808  91.340123  91.799300 109.472825 
#> 
#> $L.15
#>    4D.m.6    4D.m.5    4D.m.4    4D.m.3    4D.m.2    4D.m.1 
#>  0.000000  4.600091  8.776795 10.651111 12.143891 13.255588 
#> 
#> $L.16
#>     5A.m.54     5A.m.53     5A.m.52     5A.m.51     5A.m.50     5A.m.49 
#>   0.0000000   0.4587285   0.9174569   1.8376979  19.5108820  38.7568771 
#>     5A.m.48     5A.m.47     5A.m.46     5A.m.45     5A.m.44     5A.m.43 
#>  39.2156055  41.0512923  49.8568071  50.3155356  50.7742640  51.2329925 
#>     5A.m.42     5A.m.41     5A.m.40     5A.m.39     5A.m.38     5A.m.36 
#>  51.6917210  53.0682154  58.1313220  58.5900505  60.4257372  61.8022316 
#>     5A.m.37     5A.m.35     5A.m.34     5A.m.33     5A.m.32     5A.m.31 
#>  61.8022316  62.2609601  63.1784942  66.3939287  69.1490072  70.0665414 
#>     5A.m.30     5A.m.29     5A.m.28     5A.m.27     5A.m.26     5A.m.25 
#>  70.5252698  72.3609565  73.7374510  74.1961794  76.0318662  79.7082016 
#>     5A.m.24     5A.m.23     5A.m.22     5A.m.21     5A.m.20     5A.m.19 
#>  80.1669301  86.1587514  87.0762856 100.7069618 105.7700684 106.2287969 
#>     5A.m.18     5A.m.17     5A.m.16     5A.m.15     5A.m.14     5A.m.13 
#> 106.6875254 108.0640198 109.4405142 117.7737037 119.6093904 123.7472514 
#>     5A.m.12     5A.m.10     5A.m.11      5A.m.9      5A.m.8      5A.m.7 
#> 128.3473425 130.6425312 131.1012597 160.1198614 166.1116828 171.1747894 
#>      5A.m.6      5A.m.5      5A.m.4      5A.m.3      5A.m.2      5A.m.1 
#> 172.0923235 172.5510520 173.0097805 174.8454672 175.3041956 181.7618822 
#> 
#> $L.17
#>     5B.m.56     5B.m.55     5B.m.54     5B.m.53     5B.m.52     5B.m.51 
#>   0.0000000   0.9175342  14.0558285  17.7321639  20.0273526  20.4860811 
#>     5B.m.50     5B.m.49     5B.m.48     5B.m.47     5B.m.46     5B.m.45 
#>  22.3217678  22.7804963  23.2392247  23.6979532  24.1566817  24.6154101 
#>     5B.m.44     5B.m.43     5B.m.42     5B.m.41     5B.m.40     5B.m.39 
#>  25.0741386  25.9916728  26.4504012  26.9091297  27.3678582  27.8265933 
#>     5B.m.38     5B.m.37     5B.m.36     5B.m.35     5B.m.34     5B.m.33 
#>  28.7680562  29.2508125  31.5461135  32.0048420  32.4635705  33.3811046 
#>     5B.m.32     5B.m.31     5B.m.29     5B.m.30     5B.m.28     5B.m.27 
#>  33.8398331  34.7573673  36.1338617  36.5925902  37.0513186  60.1193592 
#>     5B.m.26     5B.m.25     5B.m.24     5B.m.23     5B.m.22     5B.m.21 
#>  64.2572202  64.7159486  65.1746771  66.5511715  72.0781607  75.2935951 
#>     5B.m.20     5B.m.19     5B.m.18     5B.m.17     5B.m.15     5B.m.16 
#>  76.7000704  78.6053647  79.5624655  83.2388010  83.6975294  83.6975294 
#>     5B.m.14     5B.m.13     5B.m.12     5B.m.11     5B.m.10      5B.m.9 
#>  98.3208693  98.7795978 111.9178921 112.8354262 115.5905047 120.6536113 
#>      5B.m.8      5B.m.7      5B.m.6      5B.m.5      5B.m.4      5B.m.3 
#> 125.2537024 125.7124309 135.4673718 137.7625605 138.6800947 139.6333935 
#>      5B.m.2      5B.m.1 
#> 142.1178003 142.7883741 
#> 
#> $L.18
#>     5D.m.6     5D.m.4     5D.m.5     5D.m.1     5D.m.2     5D.m.3 
#>  0.0000000  0.4587285  0.4587285 17.1024336 17.5611621 18.0198906 
#> 
#> $L.19
#>    6A.m.27    6A.m.26    6A.m.25    6A.m.24    6A.m.23    6A.m.22    6A.m.21 
#>  0.0000000  0.4587285  3.2138070  3.6725354  4.5900696 29.3742386 37.2365718 
#>    6A.m.20    6A.m.18    6A.m.19    6A.m.17    6A.m.16    6A.m.15    6A.m.14 
#> 38.6130662 40.9082550 41.3669834 41.8257119 42.2844404 43.2019745 45.0376613 
#>    6A.m.13    6A.m.12    6A.m.11    6A.m.10     6A.m.9     6A.m.8     6A.m.7 
#> 51.9623308 56.1001918 63.9625251 65.3494460 74.1657062 74.6244347 79.6875413 
#>     6A.m.6     6A.m.5     6A.m.4     6A.m.3     6A.m.2     6A.m.1 
#> 91.3631163 93.1988030 95.0344897 97.3296785 98.2472126 98.7059411 
#> 
#> $L.2
#>  1B1.m.4  1B1.m.5  1B1.m.3  1B1.m.2  1B1.m.1 
#> 0.000000 0.000000 4.137861 6.892939 7.810474 
#> 
#> $L.20
#>    6B.m.41    6B.m.40    6B.m.39    6B.m.38    6B.m.37    6B.m.36    6B.m.35 
#>   0.000000   3.228588   4.666779   7.462087   9.716855  10.674985  11.133714 
#>    6B.m.34    6B.m.33    6B.m.32    6B.m.31    6B.m.30    6B.m.29    6B.m.28 
#>  13.428903  20.821759  23.152283  26.403125  43.559982  46.315061  46.773789 
#>    6B.m.27    6B.m.26    6B.m.25    6B.m.24    6B.m.23    6B.m.22    6B.m.21 
#>  47.232518  47.694728  55.553493  56.012221  58.767300  60.143794  60.602523 
#>    6B.m.20    6B.m.19    6B.m.18    6B.m.16    6B.m.17    6B.m.15    6B.m.14 
#>  62.438210  63.923978  64.732295  65.649829  66.108558  66.567286  67.943781 
#>    6B.m.13    6B.m.12    6B.m.11    6B.m.10     6B.m.9     6B.m.8     6B.m.7 
#>  68.402509  69.320043  73.457904  80.382574  80.841302  81.300031 101.078220 
#>     6B.m.6     6B.m.5     6B.m.4     6B.m.2     6B.m.3     6B.m.1 
#> 110.357621 110.816350 112.192844 114.488173 114.947042 115.406325 
#> 
#> $L.21
#>    6D.m.15    6D.m.14    6D.m.13    6D.m.12    6D.m.11    6D.m.10     6D.m.9 
#>  0.0000000  0.9225051  4.1379396 10.6171243 39.6108758 57.8493328 58.7668670 
#>     6D.m.8     6D.m.7     6D.m.6     6D.m.5     6D.m.4     6D.m.2     6D.m.3 
#> 59.2262509 69.4591643 69.9178927 70.8354269 71.2941554 71.7528838 71.7528838 
#>     6D.m.1 
#> 72.2116123 
#> 
#> $L.22
#>    7A.m.32    7A.m.33    7A.m.31    7A.m.30    7A.m.29    7A.m.28    7A.m.27 
#>   0.000000   0.000000   0.684333   4.106999   5.024533   5.488570  10.094015 
#>    7A.m.26    7A.m.25    7A.m.24    7A.m.23    7A.m.22    7A.m.21    7A.m.20 
#>  11.011549  18.873882  34.501225  36.796414  37.715916  47.472900  68.874166 
#>    7A.m.19    7A.m.18    7A.m.17    7A.m.16    7A.m.15    7A.m.14    7A.m.13 
#>  69.791700  70.709234  71.167963  72.085497  72.544226  73.920720  78.474431 
#>    7A.m.12    7A.m.11    7A.m.10     7A.m.9     7A.m.8     7A.m.7     7A.m.6 
#>  79.437975  80.814469  83.109658  84.027192  86.782271  87.240999 113.797115 
#>     7A.m.5     7A.m.4     7A.m.3     7A.m.2     7A.m.1 
#> 118.860222 127.193412 131.793503 133.796103 134.474614 
#> 
#> $L.23
#>   7B.m.37   7B.m.36   7B.m.35   7B.m.34   7B.m.33   7B.m.31   7B.m.32   7B.m.30 
#>   0.00000  23.63430  24.09303  27.76936  29.14586  30.06339  30.52212  36.04911 
#>   7B.m.29   7B.m.28   7B.m.27   7B.m.26   7B.m.25   7B.m.24   7B.m.23   7B.m.21 
#>  38.80419  39.72172  58.43965  59.35718  62.57262  65.32769  65.78642  69.00186 
#>   7B.m.22   7B.m.20   7B.m.19   7B.m.17   7B.m.18   7B.m.16   7B.m.15   7B.m.14 
#>  69.00186  69.91939  72.67447  74.96966  75.42839  75.88712  76.80465  77.26338 
#>   7B.m.13   7B.m.12   7B.m.11   7B.m.10    7B.m.9    7B.m.8    7B.m.6    7B.m.7 
#>  79.55857  80.47610  81.39364  81.85236  82.31109  83.22863  84.60512  85.36757 
#>    7B.m.5    7B.m.4    7B.m.3    7B.m.2    7B.m.1 
#>  86.91130  92.42609  92.88482  98.87664 108.15604 
#> 
#> $L.24
#>    7D.m.4    7D.m.2    7D.m.3    7D.m.1 
#>  0.000000  5.526989  5.526989 21.660902 
#> 
#> $L.3
#>  1B2.m.33  1B2.m.32  1B2.m.31  1B2.m.30  1B2.m.29  1B2.m.28  1B2.m.27  1B2.m.26 
#>  0.000000  9.279401 10.655895 11.573429 12.032158 13.867845 15.289476 16.711107 
#>  1B2.m.25  1B2.m.24  1B2.m.23  1B2.m.22  1B2.m.21  1B2.m.20  1B2.m.19  1B2.m.18 
#> 17.169836 18.546330 19.463864 20.851965 29.173232 31.008919 32.844606 36.520941 
#>  1B2.m.17  1B2.m.16  1B2.m.15  1B2.m.14  1B2.m.13  1B2.m.12  1B2.m.11  1B2.m.10 
#> 36.979670 38.356164 38.814893 40.650579 42.027074 42.944608 43.862142 44.320871 
#>   1B2.m.9   1B2.m.8   1B2.m.7   1B2.m.6   1B2.m.5   1B2.m.4   1B2.m.3   1B2.m.2 
#> 44.779599 46.156094 47.991780 48.450509 49.827003 51.203498 51.662226 59.995416 
#>   1B2.m.1 
#> 60.912950 
#> 
#> $L.4
#>  1D.m.10   1D.m.9   1D.m.8   1D.m.7   1D.m.6   1D.m.5   1D.m.3   1D.m.4 
#>  0.00000 22.50673 22.96620 23.88404 55.53218 59.75146 78.46976 78.46976 
#>   1D.m.2   1D.m.1 
#> 78.92881 81.22400 
#> 
#> $L.5
#>    2A.m.40    2A.m.39    2A.m.38    2A.m.37    2A.m.36    2A.m.35    2A.m.34 
#>   0.000000   2.654039   3.575343   4.953997   5.872637   7.249131  18.886363 
#>    2A.m.33    2A.m.32    2A.m.31    2A.m.30    2A.m.29    2A.m.28    2A.m.26 
#>  24.338629  24.797357  29.860464  30.319193  45.946535  46.885701  49.640780 
#>    2A.m.27    2A.m.25    2A.m.24    2A.m.23    2A.m.22    2A.m.21    2A.m.20 
#>  51.476466  52.394000  54.689189  57.904624  58.363352  58.822081  62.959942 
#>    2A.m.19    2A.m.18    2A.m.17    2A.m.16    2A.m.15    2A.m.14    2A.m.13 
#>  63.418670  76.067061  77.902748  80.657826  83.412910  83.886982  84.345711 
#>    2A.m.12    2A.m.11    2A.m.10     2A.m.9     2A.m.8     2A.m.7     2A.m.6 
#>  84.804439  85.263167  85.721896 101.855809 105.993670 106.452398 111.979387 
#>     2A.m.5     2A.m.4     2A.m.3     2A.m.2     2A.m.1 
#> 126.105039 127.940726 128.858260 138.613201 153.236541 
#> 
#> $L.6
#>  2B.m.35  2B.m.33  2B.m.34  2B.m.31  2B.m.32  2B.m.30  2B.m.29  2B.m.28 
#>  0.00000 16.64371 17.10243 18.47893 18.47893 18.93766 19.39638 19.85511 
#>  2B.m.27  2B.m.26  2B.m.25  2B.m.24  2B.m.23  2B.m.22  2B.m.21  2B.m.20 
#> 20.77265 21.69018 22.14891 23.98460 24.44333 25.81982 27.19631 27.65504 
#>  2B.m.19  2B.m.18  2B.m.17  2B.m.16  2B.m.15  2B.m.14  2B.m.13  2B.m.12 
#> 29.03154 29.94907 30.86661 31.32533 31.78406 32.24279 33.61929 35.45497 
#>  2B.m.11  2B.m.10   2B.m.9   2B.m.8   2B.m.7   2B.m.6   2B.m.5   2B.m.4 
#> 36.37251 45.17802 45.63675 46.09548 47.93116 56.26435 65.54376 71.07074 
#>   2B.m.3   2B.m.2   2B.m.1 
#> 71.52947 94.59751 97.37680 
#> 
#> $L.7
#>  2D1.m.8  2D1.m.6  2D1.m.7  2D1.m.5  2D1.m.4  2D1.m.3  2D1.m.2  2D1.m.1 
#>  0.00000 13.13829 13.13829 16.77185 17.73195 18.20947 18.71471 19.20363 
#> 
#> $L.8
#>  2D2.m.13  2D2.m.12  2D2.m.11   2D2.m.9  2D2.m.10   2D2.m.8   2D2.m.7   2D2.m.6 
#>  0.000000  3.215434 10.137247 11.057585 11.516314 42.462800 46.604503 47.066011 
#>   2D2.m.5   2D2.m.4   2D2.m.3   2D2.m.2   2D2.m.1 
#> 49.363984 51.660473 59.519642 59.982755 60.441483 
#> 
#> $L.9
#>     3A.m.37     3A.m.35     3A.m.36     3A.m.34     3A.m.33     3A.m.32 
#>   0.0000000   0.4587285   0.4587285   1.8352229   2.2939514   4.5891402 
#>     3A.m.31     3A.m.30     3A.m.29     3A.m.28     3A.m.27     3A.m.26 
#>   5.5066743   7.4768940  15.9172288  26.6285944  38.3224121  44.3317651 
#>     3A.m.25     3A.m.24     3A.m.23     3A.m.22     3A.m.21     3A.m.20 
#>  50.7894517  56.7812731  57.2400015  63.6976881  64.6152223  65.5327565 
#>     3A.m.19     3A.m.16     3A.m.17     3A.m.18     3A.m.15     3A.m.14 
#>  68.7481909  82.8738429  83.7913771  84.2501056  85.6266000  86.0853285 
#>     3A.m.13     3A.m.10      3A.m.9     3A.m.11     3A.m.12      3A.m.8 
#>  86.5440569  90.6819179  92.0584123  93.4349067  95.7300955 126.0286198 
#>      3A.m.7      3A.m.6      3A.m.5      3A.m.4      3A.m.3      3A.m.2 
#> 126.9461540 133.4038406 134.3213748 137.9977120 144.9223833 145.3811118 
#>      3A.m.1 
#> 145.8398402 
#> 

## one linkage group at a time (possibly break established linkage
## groups)

test2 <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", trace = FALSE)
#> Number of linkage groups: 1
#> The size of the linkage groups are: 41   
#> The number of bins in each linkage group: 41 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 5    
#> The number of bins in each linkage group: 4  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 33   
#> The number of bins in each linkage group: 33 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 10   
#> The number of bins in each linkage group: 9  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 40   
#> The number of bins in each linkage group: 40 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 35   
#> The number of bins in each linkage group: 34 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 8    
#> The number of bins in each linkage group: 7  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 13   
#> The number of bins in each linkage group: 13 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 37   
#> The number of bins in each linkage group: 36 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 30   
#> The number of bins in each linkage group: 29 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 5  
#> Number of linkage groups: 2
#> The size of the linkage groups are: 9    21  
#> The number of bins in each linkage group: 8  21  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 32   
#> The number of bins in each linkage group: 32 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 6  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 54   
#> The number of bins in each linkage group: 53 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 56   
#> The number of bins in each linkage group: 55 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 5  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 27   
#> The number of bins in each linkage group: 27 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 41   
#> The number of bins in each linkage group: 41 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 15   
#> The number of bins in each linkage group: 14 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 33   
#> The number of bins in each linkage group: 32 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 37   
#> The number of bins in each linkage group: 36 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 4    
#> The number of bins in each linkage group: 3  
pull.map(test2)
#> $`1A`
#>     1A.m.41     1A.m.40     1A.m.39     1A.m.38     1A.m.37     1A.m.36 
#>   0.0000000   0.4587285   0.9174569   3.2126457  10.1373152  11.5138097 
#>     1A.m.35     1A.m.34     1A.m.33     1A.m.32     1A.m.31     1A.m.30 
#>  16.1139008  23.9762340  24.4349625  24.8936909  25.9826349  27.1876743 
#>     1A.m.29     1A.m.28     1A.m.27     1A.m.26     1A.m.25     1A.m.24 
#>  29.0233610  32.2387955  32.6975239  34.5332107  63.5518124  64.9283068 
#>     1A.m.23     1A.m.22     1A.m.21     1A.m.20     1A.m.19     1A.m.18 
#>  68.1437413  68.6024697  69.0611982  69.9787324  71.3552268  71.8139553 
#>     1A.m.17     1A.m.16     1A.m.15     1A.m.13     1A.m.14     1A.m.12 
#>  72.7314894  73.1902179  74.5667123  75.5029664  75.9804099  76.4391384 
#>     1A.m.11     1A.m.10      1A.m.9      1A.m.8      1A.m.7      1A.m.6 
#>  77.8156328  78.7331670  81.4882455  81.9469740  82.4057024  82.8644309 
#>      1A.m.5      1A.m.4      1A.m.3      1A.m.2      1A.m.1 
#>  83.3231594  84.6996538  99.8235155 102.5785940 103.9550885 
#> 
#> $`1B1`
#>  1B1.m.4  1B1.m.5  1B1.m.3  1B1.m.2  1B1.m.1 
#> 0.000000 0.000000 4.137861 6.892939 7.810474 
#> 
#> $`1B2`
#>  1B2.m.33  1B2.m.32  1B2.m.31  1B2.m.30  1B2.m.29  1B2.m.28  1B2.m.27  1B2.m.26 
#>  0.000000  9.279401 10.655895 11.573429 12.032158 13.867845 15.289476 16.711107 
#>  1B2.m.25  1B2.m.24  1B2.m.23  1B2.m.22  1B2.m.21  1B2.m.20  1B2.m.19  1B2.m.18 
#> 17.169836 18.546330 19.463864 20.851965 29.173232 31.008919 32.844606 36.520941 
#>  1B2.m.17  1B2.m.16  1B2.m.15  1B2.m.14  1B2.m.13  1B2.m.12  1B2.m.11  1B2.m.10 
#> 36.979670 38.356164 38.814893 40.650579 42.027074 42.944608 43.862142 44.320871 
#>   1B2.m.9   1B2.m.8   1B2.m.7   1B2.m.6   1B2.m.5   1B2.m.4   1B2.m.3   1B2.m.2 
#> 44.779599 46.156094 47.991780 48.450509 49.827003 51.203498 51.662226 59.995416 
#>   1B2.m.1 
#> 60.912950 
#> 
#> $`1D`
#>  1D.m.10   1D.m.9   1D.m.8   1D.m.7   1D.m.6   1D.m.5   1D.m.3   1D.m.4 
#>  0.00000 22.50673 22.96620 23.88404 55.53218 59.75146 78.46976 78.46976 
#>   1D.m.2   1D.m.1 
#> 78.92881 81.22400 
#> 
#> $`2A`
#>    2A.m.40    2A.m.39    2A.m.38    2A.m.37    2A.m.36    2A.m.35    2A.m.34 
#>   0.000000   2.654039   3.575343   4.953997   5.872637   7.249131  18.886363 
#>    2A.m.33    2A.m.32    2A.m.31    2A.m.30    2A.m.29    2A.m.28    2A.m.26 
#>  24.338629  24.797357  29.860464  30.319193  45.946535  46.885701  49.640780 
#>    2A.m.27    2A.m.25    2A.m.24    2A.m.23    2A.m.22    2A.m.21    2A.m.20 
#>  51.476466  52.394000  54.689189  57.904624  58.363352  58.822081  62.959942 
#>    2A.m.19    2A.m.18    2A.m.17    2A.m.16    2A.m.15    2A.m.14    2A.m.13 
#>  63.418670  76.067061  77.902748  80.657826  83.412910  83.886982  84.345711 
#>    2A.m.12    2A.m.11    2A.m.10     2A.m.9     2A.m.8     2A.m.7     2A.m.6 
#>  84.804439  85.263167  85.721896 101.855809 105.993670 106.452398 111.979387 
#>     2A.m.5     2A.m.4     2A.m.3     2A.m.2     2A.m.1 
#> 126.105039 127.940726 128.858260 138.613201 153.236541 
#> 
#> $`2B`
#>  2B.m.35  2B.m.33  2B.m.34  2B.m.31  2B.m.32  2B.m.30  2B.m.29  2B.m.28 
#>  0.00000 16.64371 17.10243 18.47893 18.47893 18.93766 19.39638 19.85511 
#>  2B.m.27  2B.m.26  2B.m.25  2B.m.24  2B.m.23  2B.m.22  2B.m.21  2B.m.20 
#> 20.77265 21.69018 22.14891 23.98460 24.44333 25.81982 27.19631 27.65504 
#>  2B.m.19  2B.m.18  2B.m.17  2B.m.16  2B.m.15  2B.m.14  2B.m.13  2B.m.12 
#> 29.03154 29.94907 30.86661 31.32533 31.78406 32.24279 33.61929 35.45497 
#>  2B.m.11  2B.m.10   2B.m.9   2B.m.8   2B.m.7   2B.m.6   2B.m.5   2B.m.4 
#> 36.37251 45.17802 45.63675 46.09548 47.93116 56.26435 65.54376 71.07074 
#>   2B.m.3   2B.m.2   2B.m.1 
#> 71.52947 94.59751 97.37680 
#> 
#> $`2D1`
#>  2D1.m.8  2D1.m.6  2D1.m.7  2D1.m.5  2D1.m.4  2D1.m.3  2D1.m.2  2D1.m.1 
#>  0.00000 13.13829 13.13829 16.77185 17.73195 18.20947 18.71471 19.20363 
#> 
#> $`2D2`
#>  2D2.m.13  2D2.m.12  2D2.m.11   2D2.m.9  2D2.m.10   2D2.m.8   2D2.m.7   2D2.m.6 
#>  0.000000  3.215434 10.137247 11.057585 11.516314 42.462800 46.604503 47.066011 
#>   2D2.m.5   2D2.m.4   2D2.m.3   2D2.m.2   2D2.m.1 
#> 49.363984 51.660473 59.519642 59.982755 60.441483 
#> 
#> $`3A`
#>     3A.m.37     3A.m.35     3A.m.36     3A.m.34     3A.m.33     3A.m.32 
#>   0.0000000   0.4587285   0.4587285   1.8352229   2.2939514   4.5891402 
#>     3A.m.31     3A.m.30     3A.m.29     3A.m.28     3A.m.27     3A.m.26 
#>   5.5066743   7.4768940  15.9172288  26.6285944  38.3224121  44.3317651 
#>     3A.m.25     3A.m.24     3A.m.23     3A.m.22     3A.m.21     3A.m.20 
#>  50.7894517  56.7812731  57.2400015  63.6976881  64.6152223  65.5327565 
#>     3A.m.19     3A.m.16     3A.m.17     3A.m.18     3A.m.15     3A.m.14 
#>  68.7481909  82.8738429  83.7913771  84.2501056  85.6266000  86.0853285 
#>     3A.m.13     3A.m.10      3A.m.9     3A.m.11     3A.m.12      3A.m.8 
#>  86.5440569  90.6819179  92.0584123  93.4349067  95.7300955 126.0286198 
#>      3A.m.7      3A.m.6      3A.m.5      3A.m.4      3A.m.3      3A.m.2 
#> 126.9461540 133.4038406 134.3213748 137.9977120 144.9223833 145.3811118 
#>      3A.m.1 
#> 145.8398402 
#> 
#> $`3B`
#>   3B.m.30   3B.m.29   3B.m.28   3B.m.27   3B.m.26   3B.m.25   3B.m.24   3B.m.23 
#>   0.00000  10.23223  22.39308  30.25542  34.39367  34.85278  38.99064  53.45212 
#>   3B.m.22   3B.m.21   3B.m.20   3B.m.19   3B.m.17   3B.m.18   3B.m.16   3B.m.15 
#>  56.35658  61.41969  66.48280  67.40033  68.31786  68.77659  69.23532  69.69405 
#>   3B.m.14   3B.m.13   3B.m.11   3B.m.12   3B.m.10    3B.m.9    3B.m.8    3B.m.7 
#>  73.83191  74.29064  74.74937  74.74937  75.20810  97.15983  98.07737 104.62584 
#>    3B.m.6    3B.m.5    3B.m.4    3B.m.3    3B.m.2    3B.m.1 
#> 129.29707 131.59226 132.05099 132.50972 133.42725 133.88598 
#> 
#> $`3D`
#>    3D.m.6    3D.m.5    3D.m.4    3D.m.3    3D.m.1    3D.m.2 
#> 0.0000000 0.4587314 5.0589782 8.2768742 8.7379031 8.7379031 
#> 
#> $`4A.1`
#>   4A.m.9   4A.m.7   4A.m.8   4A.m.6   4A.m.5   4A.m.3   4A.m.4   4A.m.2 
#> 0.000000 1.835687 1.835687 2.294415 2.753144 3.670678 4.129406 6.424595 
#>   4A.m.1 
#> 7.801089 
#> 
#> $`4A.2`
#>  4A.m.30  4A.m.29  4A.m.28  4A.m.27  4A.m.26  4A.m.25  4A.m.24  4A.m.23 
#>  0.00000  6.92467 18.60024 26.93343 43.06735 44.44384 44.90257 45.82010 
#>  4A.m.22  4A.m.21  4A.m.20  4A.m.19  4A.m.18  4A.m.17  4A.m.16  4A.m.15 
#> 46.73764 47.19637 48.11390 51.79024 53.16673 53.62546 83.92398 84.38271 
#>  4A.m.14  4A.m.13  4A.m.12  4A.m.11  4A.m.10 
#> 86.67790 94.54023 94.99896 95.91650 98.21169 
#> 
#> $`4B`
#>    4B.m.32    4B.m.31    4B.m.30    4B.m.29    4B.m.28    4B.m.27    4B.m.26 
#>   0.000000   8.333189  21.963866  22.881400  24.257894  29.784884  37.647217 
#>    4B.m.25    4B.m.24    4B.m.23    4B.m.22    4B.m.21    4B.m.20    4B.m.19 
#>  41.785078  42.243806  43.161340  50.086010  51.003544  52.380039  52.838767 
#>    4B.m.18    4B.m.17    4B.m.16    4B.m.15    4B.m.14    4B.m.13    4B.m.12 
#>  53.756301  56.051490  56.510218  58.345905  58.804634  59.263362  59.722091 
#>    4B.m.11    4B.m.10     4B.m.9     4B.m.8     4B.m.7     4B.m.6     4B.m.5 
#>  60.180819  60.639547  61.557082  62.015810  64.310999  78.436651  81.191729 
#>     4B.m.4     4B.m.3     4B.m.2     4B.m.1 
#>  83.946808  91.340123  91.799300 109.472825 
#> 
#> $`4D`
#>    4D.m.6    4D.m.5    4D.m.4    4D.m.3    4D.m.2    4D.m.1 
#>  0.000000  4.600091  8.776795 10.651111 12.143891 13.255588 
#> 
#> $`5A`
#>     5A.m.54     5A.m.53     5A.m.52     5A.m.51     5A.m.50     5A.m.49 
#>   0.0000000   0.4587285   0.9174569   1.8376979  19.5108820  38.7568771 
#>     5A.m.48     5A.m.47     5A.m.46     5A.m.45     5A.m.44     5A.m.43 
#>  39.2156055  41.0512923  49.8568071  50.3155356  50.7742640  51.2329925 
#>     5A.m.42     5A.m.41     5A.m.40     5A.m.39     5A.m.38     5A.m.36 
#>  51.6917210  53.0682154  58.1313220  58.5900505  60.4257372  61.8022316 
#>     5A.m.37     5A.m.35     5A.m.34     5A.m.33     5A.m.32     5A.m.31 
#>  61.8022316  62.2609601  63.1784942  66.3939287  69.1490072  70.0665414 
#>     5A.m.30     5A.m.29     5A.m.28     5A.m.27     5A.m.26     5A.m.25 
#>  70.5252698  72.3609565  73.7374510  74.1961794  76.0318662  79.7082016 
#>     5A.m.24     5A.m.23     5A.m.22     5A.m.21     5A.m.20     5A.m.19 
#>  80.1669301  86.1587514  87.0762856 100.7069618 105.7700684 106.2287969 
#>     5A.m.18     5A.m.17     5A.m.16     5A.m.15     5A.m.14     5A.m.13 
#> 106.6875254 108.0640198 109.4405142 117.7737037 119.6093904 123.7472514 
#>     5A.m.12     5A.m.10     5A.m.11      5A.m.9      5A.m.8      5A.m.7 
#> 128.3473425 130.6425312 131.1012597 160.1198614 166.1116828 171.1747894 
#>      5A.m.6      5A.m.5      5A.m.4      5A.m.3      5A.m.2      5A.m.1 
#> 172.0923235 172.5510520 173.0097805 174.8454672 175.3041956 181.7618822 
#> 
#> $`5B`
#>     5B.m.56     5B.m.55     5B.m.54     5B.m.53     5B.m.52     5B.m.51 
#>   0.0000000   0.9175342  14.0558285  17.7321639  20.0273526  20.4860811 
#>     5B.m.50     5B.m.49     5B.m.48     5B.m.47     5B.m.46     5B.m.45 
#>  22.3217678  22.7804963  23.2392247  23.6979532  24.1566817  24.6154101 
#>     5B.m.44     5B.m.43     5B.m.42     5B.m.41     5B.m.40     5B.m.39 
#>  25.0741386  25.9916728  26.4504012  26.9091297  27.3678582  27.8265933 
#>     5B.m.38     5B.m.37     5B.m.36     5B.m.35     5B.m.34     5B.m.33 
#>  28.7680562  29.2508125  31.5461135  32.0048420  32.4635705  33.3811046 
#>     5B.m.32     5B.m.31     5B.m.29     5B.m.30     5B.m.28     5B.m.27 
#>  33.8398331  34.7573673  36.1338617  36.5925902  37.0513186  60.1193592 
#>     5B.m.26     5B.m.25     5B.m.24     5B.m.23     5B.m.22     5B.m.21 
#>  64.2572202  64.7159486  65.1746771  66.5511715  72.0781607  75.2935951 
#>     5B.m.20     5B.m.19     5B.m.18     5B.m.17     5B.m.15     5B.m.16 
#>  76.7000704  78.6053647  79.5624655  83.2388010  83.6975294  83.6975294 
#>     5B.m.14     5B.m.13     5B.m.12     5B.m.11     5B.m.10      5B.m.9 
#>  98.3208693  98.7795978 111.9178921 112.8354262 115.5905047 120.6536113 
#>      5B.m.8      5B.m.7      5B.m.6      5B.m.5      5B.m.4      5B.m.3 
#> 125.2537024 125.7124309 135.4673718 137.7625605 138.6800947 139.6333935 
#>      5B.m.2      5B.m.1 
#> 142.1178003 142.7883741 
#> 
#> $`5D`
#>     5D.m.6     5D.m.4     5D.m.5     5D.m.1     5D.m.2     5D.m.3 
#>  0.0000000  0.4587285  0.4587285 17.1024336 17.5611621 18.0198906 
#> 
#> $`6A`
#>    6A.m.27    6A.m.26    6A.m.25    6A.m.24    6A.m.23    6A.m.22    6A.m.21 
#>  0.0000000  0.4587285  3.2138070  3.6725354  4.5900696 29.3742386 37.2365718 
#>    6A.m.20    6A.m.18    6A.m.19    6A.m.17    6A.m.16    6A.m.15    6A.m.14 
#> 38.6130662 40.9082550 41.3669834 41.8257119 42.2844404 43.2019745 45.0376613 
#>    6A.m.13    6A.m.12    6A.m.11    6A.m.10     6A.m.9     6A.m.8     6A.m.7 
#> 51.9623308 56.1001918 63.9625251 65.3494460 74.1657062 74.6244347 79.6875413 
#>     6A.m.6     6A.m.5     6A.m.4     6A.m.3     6A.m.2     6A.m.1 
#> 91.3631163 93.1988030 95.0344897 97.3296785 98.2472126 98.7059411 
#> 
#> $`6B`
#>    6B.m.41    6B.m.40    6B.m.39    6B.m.38    6B.m.37    6B.m.36    6B.m.35 
#>   0.000000   3.228588   4.666779   7.462087   9.716855  10.674985  11.133714 
#>    6B.m.34    6B.m.33    6B.m.32    6B.m.31    6B.m.30    6B.m.29    6B.m.28 
#>  13.428903  20.821759  23.152283  26.403125  43.559982  46.315061  46.773789 
#>    6B.m.27    6B.m.26    6B.m.25    6B.m.24    6B.m.23    6B.m.22    6B.m.21 
#>  47.232518  47.694728  55.553493  56.012221  58.767300  60.143794  60.602523 
#>    6B.m.20    6B.m.19    6B.m.18    6B.m.16    6B.m.17    6B.m.15    6B.m.14 
#>  62.438210  63.923978  64.732295  65.649829  66.108558  66.567286  67.943781 
#>    6B.m.13    6B.m.12    6B.m.11    6B.m.10     6B.m.9     6B.m.8     6B.m.7 
#>  68.402509  69.320043  73.457904  80.382574  80.841302  81.300031 101.078220 
#>     6B.m.6     6B.m.5     6B.m.4     6B.m.2     6B.m.3     6B.m.1 
#> 110.357621 110.816350 112.192844 114.488173 114.947042 115.406325 
#> 
#> $`6D`
#>    6D.m.15    6D.m.14    6D.m.13    6D.m.12    6D.m.11    6D.m.10     6D.m.9 
#>  0.0000000  0.9225051  4.1379396 10.6171243 39.6108758 57.8493328 58.7668670 
#>     6D.m.8     6D.m.7     6D.m.6     6D.m.5     6D.m.4     6D.m.2     6D.m.3 
#> 59.2262509 69.4591643 69.9178927 70.8354269 71.2941554 71.7528838 71.7528838 
#>     6D.m.1 
#> 72.2116123 
#> 
#> $`7A`
#>    7A.m.32    7A.m.33    7A.m.31    7A.m.30    7A.m.29    7A.m.28    7A.m.27 
#>   0.000000   0.000000   0.684333   4.106999   5.024533   5.488570  10.094015 
#>    7A.m.26    7A.m.25    7A.m.24    7A.m.23    7A.m.22    7A.m.21    7A.m.20 
#>  11.011549  18.873882  34.501225  36.796414  37.715916  47.472900  68.874166 
#>    7A.m.19    7A.m.18    7A.m.17    7A.m.16    7A.m.15    7A.m.14    7A.m.13 
#>  69.791700  70.709234  71.167963  72.085497  72.544226  73.920720  78.474431 
#>    7A.m.12    7A.m.11    7A.m.10     7A.m.9     7A.m.8     7A.m.7     7A.m.6 
#>  79.437975  80.814469  83.109658  84.027192  86.782271  87.240999 113.797115 
#>     7A.m.5     7A.m.4     7A.m.3     7A.m.2     7A.m.1 
#> 118.860222 127.193412 131.793503 133.796103 134.474614 
#> 
#> $`7B`
#>   7B.m.37   7B.m.36   7B.m.35   7B.m.34   7B.m.33   7B.m.31   7B.m.32   7B.m.30 
#>   0.00000  23.63430  24.09303  27.76936  29.14586  30.06339  30.52212  36.04911 
#>   7B.m.29   7B.m.28   7B.m.27   7B.m.26   7B.m.25   7B.m.24   7B.m.23   7B.m.21 
#>  38.80419  39.72172  58.43965  59.35718  62.57262  65.32769  65.78642  69.00186 
#>   7B.m.22   7B.m.20   7B.m.19   7B.m.17   7B.m.18   7B.m.16   7B.m.15   7B.m.14 
#>  69.00186  69.91939  72.67447  74.96966  75.42839  75.88712  76.80465  77.26338 
#>   7B.m.13   7B.m.12   7B.m.11   7B.m.10    7B.m.9    7B.m.8    7B.m.6    7B.m.7 
#>  79.55857  80.47610  81.39364  81.85236  82.31109  83.22863  84.60512  85.36757 
#>    7B.m.5    7B.m.4    7B.m.3    7B.m.2    7B.m.1 
#>  86.91130  92.42609  92.88482  98.87664 108.15604 
#> 
#> $`7D`
#>    7D.m.4    7D.m.2    7D.m.3    7D.m.1 
#>  0.000000  5.526989  5.526989 21.660902 
#> 

## one linkage group at a time (do not break established linkage groups)

test3 <- mstmap(mapDH, bychr = TRUE, dist.fun = "kosambi", p.value = 2,
                trace = FALSE)
#> Number of linkage groups: 1
#> The size of the linkage groups are: 41   
#> The number of bins in each linkage group: 41 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 5    
#> The number of bins in each linkage group: 4  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 33   
#> The number of bins in each linkage group: 33 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 10   
#> The number of bins in each linkage group: 9  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 40   
#> The number of bins in each linkage group: 40 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 35   
#> The number of bins in each linkage group: 34 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 8    
#> The number of bins in each linkage group: 7  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 13   
#> The number of bins in each linkage group: 13 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 37   
#> The number of bins in each linkage group: 36 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 30   
#> The number of bins in each linkage group: 29 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 5  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 30   
#> The number of bins in each linkage group: 29 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 32   
#> The number of bins in each linkage group: 32 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 6  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 54   
#> The number of bins in each linkage group: 53 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 56   
#> The number of bins in each linkage group: 55 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 6    
#> The number of bins in each linkage group: 5  
#> Number of linkage groups: 1
#> The size of the linkage groups are: 27   
#> The number of bins in each linkage group: 27 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 41   
#> The number of bins in each linkage group: 41 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 15   
#> The number of bins in each linkage group: 14 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 33   
#> The number of bins in each linkage group: 32 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 37   
#> The number of bins in each linkage group: 36 
#> Number of linkage groups: 1
#> The size of the linkage groups are: 4    
#> The number of bins in each linkage group: 3  
pull.map(test3)
#> $`1A`
#>     1A.m.41     1A.m.40     1A.m.39     1A.m.38     1A.m.37     1A.m.36 
#>   0.0000000   0.4587285   0.9174569   3.2126457  10.1373152  11.5138097 
#>     1A.m.35     1A.m.34     1A.m.33     1A.m.32     1A.m.31     1A.m.30 
#>  16.1139008  23.9762340  24.4349625  24.8936909  25.9826349  27.1876743 
#>     1A.m.29     1A.m.28     1A.m.27     1A.m.26     1A.m.25     1A.m.24 
#>  29.0233610  32.2387955  32.6975239  34.5332107  63.5518124  64.9283068 
#>     1A.m.23     1A.m.22     1A.m.21     1A.m.20     1A.m.19     1A.m.18 
#>  68.1437413  68.6024697  69.0611982  69.9787324  71.3552268  71.8139553 
#>     1A.m.17     1A.m.16     1A.m.15     1A.m.13     1A.m.14     1A.m.12 
#>  72.7314894  73.1902179  74.5667123  75.5029664  75.9804099  76.4391384 
#>     1A.m.11     1A.m.10      1A.m.9      1A.m.8      1A.m.7      1A.m.6 
#>  77.8156328  78.7331670  81.4882455  81.9469740  82.4057024  82.8644309 
#>      1A.m.5      1A.m.4      1A.m.3      1A.m.2      1A.m.1 
#>  83.3231594  84.6996538  99.8235155 102.5785940 103.9550885 
#> 
#> $`1B1`
#>  1B1.m.4  1B1.m.5  1B1.m.3  1B1.m.2  1B1.m.1 
#> 0.000000 0.000000 4.137861 6.892939 7.810474 
#> 
#> $`1B2`
#>  1B2.m.33  1B2.m.32  1B2.m.31  1B2.m.30  1B2.m.29  1B2.m.28  1B2.m.27  1B2.m.26 
#>  0.000000  9.279401 10.655895 11.573429 12.032158 13.867845 15.289476 16.711107 
#>  1B2.m.25  1B2.m.24  1B2.m.23  1B2.m.22  1B2.m.21  1B2.m.20  1B2.m.19  1B2.m.18 
#> 17.169836 18.546330 19.463864 20.851965 29.173232 31.008919 32.844606 36.520941 
#>  1B2.m.17  1B2.m.16  1B2.m.15  1B2.m.14  1B2.m.13  1B2.m.12  1B2.m.11  1B2.m.10 
#> 36.979670 38.356164 38.814893 40.650579 42.027074 42.944608 43.862142 44.320871 
#>   1B2.m.9   1B2.m.8   1B2.m.7   1B2.m.6   1B2.m.5   1B2.m.4   1B2.m.3   1B2.m.2 
#> 44.779599 46.156094 47.991780 48.450509 49.827003 51.203498 51.662226 59.995416 
#>   1B2.m.1 
#> 60.912950 
#> 
#> $`1D`
#>  1D.m.10   1D.m.9   1D.m.8   1D.m.7   1D.m.6   1D.m.5   1D.m.3   1D.m.4 
#>  0.00000 22.50673 22.96620 23.88404 55.53218 59.75146 78.46976 78.46976 
#>   1D.m.2   1D.m.1 
#> 78.92881 81.22400 
#> 
#> $`2A`
#>    2A.m.40    2A.m.39    2A.m.38    2A.m.37    2A.m.36    2A.m.35    2A.m.34 
#>   0.000000   2.654039   3.575343   4.953997   5.872637   7.249131  18.886363 
#>    2A.m.33    2A.m.32    2A.m.31    2A.m.30    2A.m.29    2A.m.28    2A.m.26 
#>  24.338629  24.797357  29.860464  30.319193  45.946535  46.885701  49.640780 
#>    2A.m.27    2A.m.25    2A.m.24    2A.m.23    2A.m.22    2A.m.21    2A.m.20 
#>  51.476466  52.394000  54.689189  57.904624  58.363352  58.822081  62.959942 
#>    2A.m.19    2A.m.18    2A.m.17    2A.m.16    2A.m.15    2A.m.14    2A.m.13 
#>  63.418670  76.067061  77.902748  80.657826  83.412910  83.886982  84.345711 
#>    2A.m.12    2A.m.11    2A.m.10     2A.m.9     2A.m.8     2A.m.7     2A.m.6 
#>  84.804439  85.263167  85.721896 101.855809 105.993670 106.452398 111.979387 
#>     2A.m.5     2A.m.4     2A.m.3     2A.m.2     2A.m.1 
#> 126.105039 127.940726 128.858260 138.613201 153.236541 
#> 
#> $`2B`
#>  2B.m.35  2B.m.33  2B.m.34  2B.m.31  2B.m.32  2B.m.30  2B.m.29  2B.m.28 
#>  0.00000 16.64371 17.10243 18.47893 18.47893 18.93766 19.39638 19.85511 
#>  2B.m.27  2B.m.26  2B.m.25  2B.m.24  2B.m.23  2B.m.22  2B.m.21  2B.m.20 
#> 20.77265 21.69018 22.14891 23.98460 24.44333 25.81982 27.19631 27.65504 
#>  2B.m.19  2B.m.18  2B.m.17  2B.m.16  2B.m.15  2B.m.14  2B.m.13  2B.m.12 
#> 29.03154 29.94907 30.86661 31.32533 31.78406 32.24279 33.61929 35.45497 
#>  2B.m.11  2B.m.10   2B.m.9   2B.m.8   2B.m.7   2B.m.6   2B.m.5   2B.m.4 
#> 36.37251 45.17802 45.63675 46.09548 47.93116 56.26435 65.54376 71.07074 
#>   2B.m.3   2B.m.2   2B.m.1 
#> 71.52947 94.59751 97.37680 
#> 
#> $`2D1`
#>  2D1.m.8  2D1.m.6  2D1.m.7  2D1.m.5  2D1.m.4  2D1.m.3  2D1.m.2  2D1.m.1 
#>  0.00000 13.13829 13.13829 16.77185 17.73195 18.20947 18.71471 19.20363 
#> 
#> $`2D2`
#>  2D2.m.13  2D2.m.12  2D2.m.11   2D2.m.9  2D2.m.10   2D2.m.8   2D2.m.7   2D2.m.6 
#>  0.000000  3.215434 10.137247 11.057585 11.516314 42.462800 46.604503 47.066011 
#>   2D2.m.5   2D2.m.4   2D2.m.3   2D2.m.2   2D2.m.1 
#> 49.363984 51.660473 59.519642 59.982755 60.441483 
#> 
#> $`3A`
#>     3A.m.37     3A.m.35     3A.m.36     3A.m.34     3A.m.33     3A.m.32 
#>   0.0000000   0.4587285   0.4587285   1.8352229   2.2939514   4.5891402 
#>     3A.m.31     3A.m.30     3A.m.29     3A.m.28     3A.m.27     3A.m.26 
#>   5.5066743   7.4768940  15.9172288  26.6285944  38.3224121  44.3317651 
#>     3A.m.25     3A.m.24     3A.m.23     3A.m.22     3A.m.21     3A.m.20 
#>  50.7894517  56.7812731  57.2400015  63.6976881  64.6152223  65.5327565 
#>     3A.m.19     3A.m.18     3A.m.17     3A.m.16     3A.m.15     3A.m.14 
#>  68.7481909  83.3715308  83.8302593  84.7477934  85.6653276  86.1240561 
#>     3A.m.13     3A.m.10      3A.m.9     3A.m.11     3A.m.12      3A.m.8 
#>  86.5827845  90.7206455  92.0971399  93.4736344  95.7688231 126.0673475 
#>      3A.m.7      3A.m.6      3A.m.5      3A.m.4      3A.m.3      3A.m.2 
#> 126.9848816 133.4425682 134.3601024 138.0364396 144.9611109 145.4198394 
#>      3A.m.1 
#> 145.8785678 
#> 
#> $`3B`
#>   3B.m.30   3B.m.29   3B.m.28   3B.m.27   3B.m.26   3B.m.25   3B.m.24   3B.m.23 
#>   0.00000  10.23223  22.39308  30.25542  34.39367  34.85278  38.99064  53.45212 
#>   3B.m.22   3B.m.21   3B.m.20   3B.m.19   3B.m.17   3B.m.18   3B.m.16   3B.m.15 
#>  56.35658  61.41969  66.48280  67.40033  68.31786  68.77659  69.23532  69.69405 
#>   3B.m.14   3B.m.13   3B.m.11   3B.m.12   3B.m.10    3B.m.9    3B.m.8    3B.m.7 
#>  73.83191  74.29064  74.74937  74.74937  75.20810  97.15983  98.07737 104.62584 
#>    3B.m.6    3B.m.5    3B.m.4    3B.m.3    3B.m.2    3B.m.1 
#> 129.29707 131.59226 132.05099 132.50972 133.42725 133.88598 
#> 
#> $`3D`
#>    3D.m.6    3D.m.5    3D.m.4    3D.m.3    3D.m.1    3D.m.2 
#> 0.0000000 0.4587314 5.0589782 8.2768742 8.7379031 8.7379031 
#> 
#> $`4A`
#>   4A.m.30   4A.m.29   4A.m.28   4A.m.27   4A.m.26   4A.m.25   4A.m.24   4A.m.23 
#>   0.00000   6.92467  18.60024  26.93343  43.06735  44.44384  44.90257  45.82010 
#>   4A.m.22   4A.m.21   4A.m.20   4A.m.19   4A.m.18   4A.m.17   4A.m.16   4A.m.15 
#>  46.73764  47.19637  48.11390  51.79024  53.16673  53.62546  83.92398  84.38271 
#>   4A.m.14   4A.m.13   4A.m.12   4A.m.11   4A.m.10    4A.m.9    4A.m.7    4A.m.8 
#>  86.67790  94.54023  94.99896  95.91650  98.21169 140.42607 142.26175 142.26175 
#>    4A.m.6    4A.m.5    4A.m.3    4A.m.4    4A.m.2    4A.m.1 
#> 142.72048 143.17921 144.09674 144.55547 146.85066 148.22716 
#> 
#> $`4B`
#>    4B.m.32    4B.m.31    4B.m.30    4B.m.29    4B.m.28    4B.m.27    4B.m.26 
#>   0.000000   8.333189  21.963866  22.881400  24.257894  29.784884  37.647217 
#>    4B.m.25    4B.m.24    4B.m.23    4B.m.22    4B.m.21    4B.m.20    4B.m.19 
#>  41.785078  42.243806  43.161340  50.086010  51.003544  52.380039  52.838767 
#>    4B.m.18    4B.m.17    4B.m.16    4B.m.15    4B.m.14    4B.m.13    4B.m.12 
#>  53.756301  56.051490  56.510218  58.345905  58.804634  59.263362  59.722091 
#>    4B.m.11    4B.m.10     4B.m.9     4B.m.8     4B.m.7     4B.m.6     4B.m.5 
#>  60.180819  60.639547  61.557082  62.015810  64.310999  78.436651  81.191729 
#>     4B.m.4     4B.m.3     4B.m.2     4B.m.1 
#>  83.946808  91.340123  91.799300 109.472825 
#> 
#> $`4D`
#>    4D.m.6    4D.m.5    4D.m.4    4D.m.3    4D.m.2    4D.m.1 
#>  0.000000  4.600091  8.776795 10.651111 12.143891 13.255588 
#> 
#> $`5A`
#>     5A.m.54     5A.m.53     5A.m.52     5A.m.51     5A.m.50     5A.m.49 
#>   0.0000000   0.4587285   0.9174569   1.8376979  19.5108820  38.7568771 
#>     5A.m.48     5A.m.47     5A.m.46     5A.m.45     5A.m.44     5A.m.43 
#>  39.2156055  41.0512923  49.8568071  50.3155356  50.7742640  51.2329925 
#>     5A.m.42     5A.m.41     5A.m.40     5A.m.39     5A.m.38     5A.m.36 
#>  51.6917210  53.0682154  58.1313220  58.5900505  60.4257372  61.8022316 
#>     5A.m.37     5A.m.35     5A.m.34     5A.m.33     5A.m.32     5A.m.31 
#>  61.8022316  62.2609601  63.1784942  66.3939287  69.1490072  70.0665414 
#>     5A.m.30     5A.m.29     5A.m.28     5A.m.27     5A.m.26     5A.m.25 
#>  70.5252698  72.3609565  73.7374510  74.1961794  76.0318662  79.7082016 
#>     5A.m.24     5A.m.23     5A.m.22     5A.m.21     5A.m.20     5A.m.19 
#>  80.1669301  86.1587514  87.0762856 100.7069618 105.7700684 106.2287969 
#>     5A.m.18     5A.m.17     5A.m.16     5A.m.15     5A.m.14     5A.m.13 
#> 106.6875254 108.0640198 109.4405142 117.7737037 119.6093904 123.7472514 
#>     5A.m.12     5A.m.10     5A.m.11      5A.m.9      5A.m.8      5A.m.7 
#> 128.3473425 130.6425312 131.1012597 160.1198614 166.1116828 171.1747894 
#>      5A.m.6      5A.m.5      5A.m.4      5A.m.3      5A.m.2      5A.m.1 
#> 172.0923235 172.5510520 173.0097805 174.8454672 175.3041956 181.7618822 
#> 
#> $`5B`
#>     5B.m.56     5B.m.55     5B.m.54     5B.m.53     5B.m.52     5B.m.51 
#>   0.0000000   0.9175342  14.0558285  17.7321639  20.0273526  20.4860811 
#>     5B.m.50     5B.m.49     5B.m.48     5B.m.47     5B.m.46     5B.m.45 
#>  22.3217678  22.7804963  23.2392247  23.6979532  24.1566817  24.6154101 
#>     5B.m.44     5B.m.43     5B.m.42     5B.m.41     5B.m.40     5B.m.39 
#>  25.0741386  25.9916728  26.4504012  26.9091297  27.3678582  27.8265933 
#>     5B.m.38     5B.m.37     5B.m.36     5B.m.35     5B.m.34     5B.m.33 
#>  28.7680562  29.2508125  31.5461135  32.0048420  32.4635705  33.3811046 
#>     5B.m.32     5B.m.31     5B.m.29     5B.m.30     5B.m.28     5B.m.27 
#>  33.8398331  34.7573673  36.1338617  36.5925902  37.0513186  60.1193592 
#>     5B.m.26     5B.m.25     5B.m.24     5B.m.23     5B.m.22     5B.m.21 
#>  64.2572202  64.7159486  65.1746771  66.5511715  72.0781607  75.2935951 
#>     5B.m.20     5B.m.19     5B.m.18     5B.m.17     5B.m.15     5B.m.16 
#>  76.7000704  78.6053647  79.5624655  83.2388010  83.6975294  83.6975294 
#>     5B.m.14     5B.m.13     5B.m.12     5B.m.11     5B.m.10      5B.m.9 
#>  98.3208693  98.7795978 111.9178921 112.8354262 115.5905047 120.6536113 
#>      5B.m.8      5B.m.7      5B.m.6      5B.m.5      5B.m.4      5B.m.3 
#> 125.2537024 125.7124309 135.4673718 137.7625605 138.6800947 139.6333935 
#>      5B.m.2      5B.m.1 
#> 142.1178003 142.7883741 
#> 
#> $`5D`
#>     5D.m.6     5D.m.4     5D.m.5     5D.m.1     5D.m.2     5D.m.3 
#>  0.0000000  0.4587285  0.4587285 17.1024336 17.5611621 18.0198906 
#> 
#> $`6A`
#>    6A.m.27    6A.m.26    6A.m.25    6A.m.24    6A.m.23    6A.m.22    6A.m.21 
#>  0.0000000  0.4587285  3.2138070  3.6725354  4.5900696 29.3742386 37.2365718 
#>    6A.m.20    6A.m.18    6A.m.19    6A.m.17    6A.m.16    6A.m.15    6A.m.14 
#> 38.6130662 40.9082550 41.3669834 41.8257119 42.2844404 43.2019745 45.0376613 
#>    6A.m.13    6A.m.12    6A.m.11    6A.m.10     6A.m.9     6A.m.8     6A.m.7 
#> 51.9623308 56.1001918 63.9625251 65.3494460 74.1657062 74.6244347 79.6875413 
#>     6A.m.6     6A.m.5     6A.m.4     6A.m.3     6A.m.2     6A.m.1 
#> 91.3631163 93.1988030 95.0344897 97.3296785 98.2472126 98.7059411 
#> 
#> $`6B`
#>    6B.m.41    6B.m.40    6B.m.39    6B.m.38    6B.m.37    6B.m.36    6B.m.35 
#>   0.000000   3.228588   4.666779   7.462087   9.716855  10.674985  11.133714 
#>    6B.m.34    6B.m.33    6B.m.32    6B.m.31    6B.m.30    6B.m.29    6B.m.28 
#>  13.428903  20.821759  23.152283  26.403125  43.559982  46.315061  46.773789 
#>    6B.m.27    6B.m.26    6B.m.25    6B.m.24    6B.m.23    6B.m.22    6B.m.21 
#>  47.232518  47.694728  55.553493  56.012221  58.767300  60.143794  60.602523 
#>    6B.m.20    6B.m.19    6B.m.18    6B.m.16    6B.m.17    6B.m.15    6B.m.14 
#>  62.438210  63.923978  64.732295  65.649829  66.108558  66.567286  67.943781 
#>    6B.m.13    6B.m.12    6B.m.11    6B.m.10     6B.m.9     6B.m.8     6B.m.7 
#>  68.402509  69.320043  73.457904  80.382574  80.841302  81.300031 101.078220 
#>     6B.m.6     6B.m.5     6B.m.4     6B.m.2     6B.m.3     6B.m.1 
#> 110.357621 110.816350 112.192844 114.488173 114.947042 115.406325 
#> 
#> $`6D`
#>    6D.m.15    6D.m.14    6D.m.13    6D.m.12    6D.m.11    6D.m.10     6D.m.9 
#>  0.0000000  0.9225051  4.1379396 10.6171243 39.6108758 57.8493328 58.7668670 
#>     6D.m.8     6D.m.7     6D.m.6     6D.m.5     6D.m.4     6D.m.2     6D.m.3 
#> 59.2262509 69.4591643 69.9178927 70.8354269 71.2941554 71.7528838 71.7528838 
#>     6D.m.1 
#> 72.2116123 
#> 
#> $`7A`
#>    7A.m.32    7A.m.33    7A.m.31    7A.m.30    7A.m.29    7A.m.28    7A.m.27 
#>   0.000000   0.000000   0.684333   4.106999   5.024533   5.488570  10.094015 
#>    7A.m.26    7A.m.25    7A.m.24    7A.m.23    7A.m.22    7A.m.21    7A.m.20 
#>  11.011549  18.873882  34.501225  36.796414  37.715916  47.472900  68.874166 
#>    7A.m.19    7A.m.18    7A.m.17    7A.m.16    7A.m.15    7A.m.14    7A.m.13 
#>  69.791700  70.709234  71.167963  72.085497  72.544226  73.920720  78.474431 
#>    7A.m.12    7A.m.11    7A.m.10     7A.m.9     7A.m.8     7A.m.7     7A.m.6 
#>  79.437975  80.814469  83.109658  84.027192  86.782271  87.240999 113.797115 
#>     7A.m.5     7A.m.4     7A.m.3     7A.m.2     7A.m.1 
#> 118.860222 127.193412 131.793503 133.796103 134.474614 
#> 
#> $`7B`
#>   7B.m.37   7B.m.36   7B.m.35   7B.m.34   7B.m.33   7B.m.31   7B.m.32   7B.m.30 
#>   0.00000  23.63430  24.09303  27.76936  29.14586  30.06339  30.52212  36.04911 
#>   7B.m.29   7B.m.28   7B.m.27   7B.m.26   7B.m.25   7B.m.24   7B.m.23   7B.m.21 
#>  38.80419  39.72172  58.43965  59.35718  62.57262  65.32769  65.78642  69.00186 
#>   7B.m.22   7B.m.20   7B.m.19   7B.m.17   7B.m.18   7B.m.16   7B.m.15   7B.m.14 
#>  69.00186  69.91939  72.67447  74.96966  75.42839  75.88712  76.80465  77.26338 
#>   7B.m.13   7B.m.12   7B.m.11   7B.m.10    7B.m.9    7B.m.8    7B.m.6    7B.m.7 
#>  79.55857  80.47610  81.39364  81.85236  82.31109  83.22863  84.60512  85.36757 
#>    7B.m.5    7B.m.4    7B.m.3    7B.m.2    7B.m.1 
#>  86.91130  92.42609  92.88482  98.87664 108.15604 
#> 
#> $`7D`
#>    7D.m.4    7D.m.2    7D.m.3    7D.m.1 
#>  0.000000  5.526989  5.526989 21.660902 
#> 

## impute before clustering and detect genotyping errors, pipe output to
## file

test4 <- mstmap(mapDH, bychr = FALSE, dist.fun = "kosambi",
                trace = TRUE, mvest.bc = TRUE, detectBadData = TRUE)
pull.map(test4)
#> $L.1
#>    1A.m.41    1A.m.40    1A.m.39    1A.m.38    1A.m.36    1A.m.37    1A.m.35 
#>  0.0000000  0.0000000  0.0000000  0.4707376  9.9782181 10.3684549 15.1908593 
#>    1A.m.34    1A.m.33    1A.m.32    1A.m.31    1A.m.30    1A.m.29    1A.m.27 
#> 22.5837154 22.5837154 23.0424439 23.0424439 23.0424439 23.5243019 28.5991334 
#>    1A.m.28    1A.m.26    1A.m.25    1A.m.24    1A.m.23    1A.m.22    1A.m.21 
#> 28.5991334 29.0002364 62.6659641 62.6659641 66.3422926 66.3422926 66.3422926 
#>    1A.m.20    1A.m.19    1A.m.18    1A.m.17    1A.m.16    1A.m.15    1A.m.13 
#> 67.7166334 69.5510783 69.5510783 70.4686124 70.4686124 72.7611026 72.7611026 
#>    1A.m.14    1A.m.12    1A.m.11    1A.m.10     1A.m.9     1A.m.8     1A.m.7 
#> 72.7611026 72.7611026 74.5965083 74.5965083 78.7322150 78.7322150 79.1909435 
#>     1A.m.6     1A.m.5     1A.m.4     1A.m.3     1A.m.1     1A.m.2 
#> 79.6496719 79.6496719 80.1096325 99.8815346 99.8815346 99.8815346 
#> 
#> $L.10
#>   3B.m.30   3B.m.29   3B.m.28   3B.m.27   3B.m.25   3B.m.26   3B.m.24   3B.m.22 
#>   0.00000   0.00000  19.45442  22.91261  23.38085  23.38085  23.81417  45.76099 
#>   3B.m.23   3B.m.21   3B.m.20   3B.m.19   3B.m.18   3B.m.17   3B.m.16   3B.m.15 
#>  45.76099  50.36004  56.35033  56.80465  56.80465  56.80465  57.26338  57.26338 
#>   3B.m.14   3B.m.11   3B.m.12   3B.m.10   3B.m.13    3B.m.8    3B.m.9    3B.m.7 
#>  62.32647  62.78521  62.78521  62.78521  63.24394  85.23370  85.23370  86.15508 
#>    3B.m.6    3B.m.5    3B.m.4    3B.m.3    3B.m.1    3B.m.2 
#> 121.91725 122.37604 122.37604 122.37604 123.74933 123.74933 
#> 
#> $L.11
#>   3D.m.5   3D.m.6   3D.m.4   3D.m.3   3D.m.1   3D.m.2 
#> 0.000000 0.000000 4.579348 7.794824 7.794824 7.794824 
#> 
#> $L.12
#>    4A.m.7    4A.m.8    4A.m.9    4A.m.6    4A.m.5    4A.m.4    4A.m.3    4A.m.2 
#> 0.0000000 0.0000000 0.4295604 0.4295604 0.8882888 1.3470173 1.8057458 4.4953910 
#>    4A.m.1 
#> 4.6067660 
#> 
#> $L.13
#>  4A.m.29  4A.m.30  4A.m.28  4A.m.27  4A.m.26  4A.m.25  4A.m.24  4A.m.23 
#>  0.00000  0.00000 10.74765 17.68147 33.81539 34.27024 34.27024 35.18777 
#>  4A.m.21  4A.m.22  4A.m.20  4A.m.19  4A.m.18  4A.m.17  4A.m.15  4A.m.16 
#> 36.23123 36.42670 36.73804 41.79060 41.79060 41.79060 72.78102 72.78102 
#>  4A.m.14  4A.m.12  4A.m.13  4A.m.11  4A.m.10 
#> 73.22419 83.62258 83.85181 84.19123 84.99929 
#> 
#> $L.14
#>  4B.m.32  4B.m.31  4B.m.30  4B.m.29  4B.m.28  4B.m.27  4B.m.26  4B.m.24 
#>  0.00000  0.00000 13.13829 13.13829 13.13829 19.59578 30.30334 30.47647 
#>  4B.m.25  4B.m.23  4B.m.22  4B.m.21  4B.m.19  4B.m.20  4B.m.18  4B.m.16 
#> 30.70000 30.95475 39.75684 39.75684 41.29357 41.51076 41.81285 45.09235 
#>  4B.m.17  4B.m.15  4B.m.14  4B.m.13  4B.m.12  4B.m.11  4B.m.10   4B.m.9 
#> 45.24587 47.90184 47.90184 48.36057 48.81930 49.27803 49.27803 50.64872 
#>   4B.m.8   4B.m.7   4B.m.6   4B.m.5   4B.m.4   4B.m.2   4B.m.3   4B.m.1 
#> 50.64872 51.10700 69.82403 69.82403 69.82403 80.04500 80.04500 80.50391 
#> 
#> $L.15
#>   4D.m.6   4D.m.5   4D.m.4   4D.m.3   4D.m.2   4D.m.1 
#> 0.000000 0.000000 5.989233 5.989233 7.501043 7.694335 
#> 
#> $L.16
#>    5A.m.54    5A.m.53    5A.m.52    5A.m.51    5A.m.50    5A.m.48    5A.m.49 
#>   0.000000   0.000000   0.000000   0.450987  17.094394  35.352126  35.455101 
#>    5A.m.47    5A.m.46    5A.m.45    5A.m.44    5A.m.43    5A.m.42    5A.m.41 
#>  35.868137  47.043409  47.043409  47.502137  47.960866  47.960866  48.428013 
#>    5A.m.39    5A.m.40    5A.m.38    5A.m.36    5A.m.37    5A.m.35    5A.m.34 
#>  54.761737  54.932223  58.026116  58.490930  58.490930  58.490930  58.926609 
#>    5A.m.33    5A.m.32    5A.m.31    5A.m.30    5A.m.29    5A.m.28    5A.m.27 
#>  63.062247  66.737951  66.737951  66.737951  68.995831  69.913366  69.913366 
#>    5A.m.26    5A.m.24    5A.m.25    5A.m.23    5A.m.22    5A.m.21    5A.m.20 
#>  69.913366  75.443389  75.443389  81.900131  81.900131 102.206118 102.664841 
#>    5A.m.19    5A.m.18    5A.m.17    5A.m.16    5A.m.14    5A.m.15    5A.m.13 
#> 102.664841 102.664841 104.500194 104.500194 114.778088 115.519843 120.975899 
#>    5A.m.12    5A.m.11    5A.m.10     5A.m.9     5A.m.8     5A.m.7     5A.m.6 
#> 127.430562 127.430562 127.430562 163.995004 168.813300 169.951103 170.420974 
#>     5A.m.5     5A.m.4     5A.m.3     5A.m.2     5A.m.1 
#> 170.420974 170.420974 172.716777 172.716777 172.716777 
#> 
#> $L.17
#>   5B.m.55   5B.m.56   5B.m.54   5B.m.53   5B.m.52   5B.m.51   5B.m.50   5B.m.49 
#>   0.00000   0.00000  18.15690  20.41107  20.54125  20.68951  23.34758  23.34758 
#>   5B.m.48   5B.m.47   5B.m.46   5B.m.45   5B.m.44   5B.m.43   5B.m.42   5B.m.41 
#>  23.80631  24.26504  24.72377  25.18250  25.18250  27.00921  27.00921  27.46794 
#>   5B.m.40   5B.m.39   5B.m.38   5B.m.37   5B.m.36   5B.m.35   5B.m.34   5B.m.33 
#>  27.92667  27.92667  29.30019  29.30019  32.51499  32.51499  32.97372  33.89125 
#>   5B.m.32   5B.m.31   5B.m.30   5B.m.29   5B.m.28   5B.m.27   5B.m.26   5B.m.25 
#>  33.89125  35.26078  36.17831  36.17831  36.17831  65.19352  65.65469  65.65469 
#>   5B.m.24   5B.m.23   5B.m.22   5B.m.21   5B.m.20   5B.m.19   5B.m.18   5B.m.17 
#>  65.65469  66.11708  73.03634  77.63305  77.63305  78.55295  78.55295  82.69080 
#>   5B.m.15   5B.m.16   5B.m.13   5B.m.14   5B.m.11   5B.m.12   5B.m.10    5B.m.9 
#>  82.69080  82.69080  97.31725  97.31725 111.14453 111.44923 112.75248 120.01552 
#>    5B.m.7    5B.m.8    5B.m.6    5B.m.5    5B.m.4    5B.m.3    5B.m.2    5B.m.1 
#> 124.62152 124.62152 136.27801 137.19685 137.19685 137.19685 140.61763 140.87587 
#> 
#> $L.18
#>   5D.m.6   5D.m.4   5D.m.5   5D.m.3   5D.m.2   5D.m.1 
#>  0.00000  0.00000  0.00000 17.15686 17.15686 17.15686 
#> 
#> $L.19
#>   6A.m.26   6A.m.27   6A.m.24   6A.m.25   6A.m.23   6A.m.21   6A.m.22   6A.m.20 
#>  0.000000  0.000000  3.345266  3.574243  3.805474 38.894390 38.894390 38.894390 
#>   6A.m.19   6A.m.18   6A.m.17   6A.m.16   6A.m.15   6A.m.14   6A.m.13   6A.m.12 
#> 40.730077 40.730077 41.188805 41.188805 41.749814 42.553294 51.358808 51.358808 
#>   6A.m.10   6A.m.11    6A.m.8    6A.m.9    6A.m.7    6A.m.6    6A.m.5    6A.m.4 
#> 61.546796 61.546796 68.948092 68.948092 69.453075 86.547562 86.547562 88.383248 
#>    6A.m.3    6A.m.1    6A.m.2 
#> 91.597386 91.597386 91.597386 
#> 
#> $L.2
#>  1B1.m.4  1B1.m.5  1B1.m.3  1B1.m.1  1B1.m.2 
#> 0.000000 0.000000 4.137861 6.892939 6.892939 
#> 
#> $L.20
#>     6B.m.39     6B.m.40     6B.m.41     6B.m.37     6B.m.35     6B.m.36 
#>   0.0000000   0.4515493   1.5257130   2.9038432   2.9038432   3.0068008 
#>     6B.m.34     6B.m.38     6B.m.33     6B.m.31     6B.m.32     6B.m.30 
#>   3.4044709   3.8879560  15.5711836  15.5711836  15.5711836  39.7776891 
#>     6B.m.29     6B.m.28     6B.m.27     6B.m.26     6B.m.24     6B.m.25 
#>  40.2362277  40.2362277  40.6949562  40.6949562  49.0971561  49.2125657 
#>     6B.m.23     6B.m.22     6B.m.21     6B.m.20     6B.m.17     6B.m.19 
#>  53.2816940  53.2816940  53.2816940  54.1992282  54.1992282  54.1992282 
#>     6B.m.18     6B.m.16     6B.m.15     6B.m.13     6B.m.14     6B.m.12 
#>  54.1992282  54.1992282  54.1992282  56.1629237  56.3624912  56.6866103 
#>     6B.m.11     6B.m.10      6B.m.9      6B.m.8      6B.m.7      6B.m.5 
#>  57.5766755  68.2880410  68.2880410  68.2880410 100.5596715 100.6728366 
#>      6B.m.6      6B.m.4      6B.m.2      6B.m.3      6B.m.1 
#> 100.8327267 102.5641833 104.3998701 104.3998701 104.3998701 
#> 
#> $L.21
#>    6D.m.14    6D.m.15    6D.m.13    6D.m.12    6D.m.11    6D.m.10     6D.m.9 
#>  0.0000000  0.2648916  1.2563401  3.8123288 60.9675620 62.3507222 62.3507222 
#>     6D.m.8     6D.m.7     6D.m.6     6D.m.5     6D.m.4     6D.m.1     6D.m.2 
#> 62.3507222 73.5420794 73.5420794 74.0008078 74.0008078 74.4595363 74.4595363 
#>     6D.m.3 
#> 74.4595363 
#> 
#> $L.22
#>     7A.m.32     7A.m.33     7A.m.31     7A.m.30     7A.m.29     7A.m.28 
#>   0.0000000   0.0000000   0.2304851   5.0439148   5.0439148   5.0439148 
#>     7A.m.26     7A.m.27     7A.m.25     7A.m.24     7A.m.23     7A.m.22 
#>  10.1179695  10.1179695  11.0321717  36.3762789  36.3762789  36.3762789 
#>     7A.m.21     7A.m.20     7A.m.19     7A.m.18     7A.m.17     7A.m.16 
#>  46.5955942  68.5473322  68.9993429  68.9993429  68.9993429  69.9168771 
#>     7A.m.15     7A.m.14     7A.m.13     7A.m.12     7A.m.11     7A.m.10 
#>  69.9168771  70.3797639  77.2967170  77.2967170  78.6732115  80.0497109 
#>      7A.m.9      7A.m.8      7A.m.7      7A.m.6      7A.m.5      7A.m.4 
#>  80.0497109  82.8047894  82.8047894 117.1717214 117.1717214 130.1204011 
#>      7A.m.3      7A.m.2      7A.m.1 
#> 132.1081693 132.2517429 132.4461639 
#> 
#> $L.23
#>    7B.m.36    7B.m.35    7B.m.37    7B.m.34    7B.m.33    7B.m.32    7B.m.31 
#>  0.0000000  0.0000000  0.5569769  5.5455454  6.0212678  6.0212678  6.0212678 
#>    7B.m.30    7B.m.29    7B.m.28    7B.m.26    7B.m.27    7B.m.25    7B.m.23 
#> 11.5482570 14.3033355 14.3033355 34.0815252 34.4196216 38.3284694 41.1219132 
#>    7B.m.24    7B.m.21    7B.m.22    7B.m.20    7B.m.19    7B.m.18    7B.m.17 
#> 41.1219132 44.7595183 44.7595183 44.7595183 50.2699348 50.2699348 50.2699348 
#>    7B.m.16    7B.m.15    7B.m.14    7B.m.13    7B.m.12    7B.m.11    7B.m.10 
#> 50.2699348 51.6428632 51.6428632 55.3154478 55.3154478 56.6890272 56.6890272 
#>     7B.m.9     7B.m.8     7B.m.7     7B.m.6     7B.m.5     7B.m.3     7B.m.4 
#> 56.6890272 58.0595921 58.6967364 58.9806390 59.4429459 66.3725092 66.3725092 
#>     7B.m.2     7B.m.1 
#> 66.9559701 73.2819914 
#> 
#> $L.24
#> 7D.m.4 7D.m.2 7D.m.3 7D.m.1 
#>      0      0      0      0 
#> 
#> $L.3
#>   1B2.m.33   1B2.m.32   1B2.m.31   1B2.m.30   1B2.m.29   1B2.m.28   1B2.m.27 
#>  0.0000000  0.4587285  1.3029093  1.4289748  1.4289748  4.1210023  4.1210023 
#>   1B2.m.26   1B2.m.25   1B2.m.24   1B2.m.23   1B2.m.22   1B2.m.21   1B2.m.20 
#>  4.5826840  4.5826840  6.4173981  6.4173981  7.2897289 16.0941551 16.0941551 
#>   1B2.m.19   1B2.m.17   1B2.m.18   1B2.m.16   1B2.m.15   1B2.m.14   1B2.m.13 
#> 17.9298418 21.7593855 21.9655462 24.0996427 24.0996427 24.5583711 24.5583711 
#>   1B2.m.12   1B2.m.11   1B2.m.10    1B2.m.9    1B2.m.8    1B2.m.6    1B2.m.7 
#> 25.4759053 26.8476927 26.8476927 26.8476927 28.6831147 30.6073916 30.7516149 
#>    1B2.m.5    1B2.m.4    1B2.m.3    1B2.m.1    1B2.m.2 
#> 32.4964175 33.8729120 33.8729120 42.6784268 42.6784268 
#> 
#> $L.4
#>     1D.m.8     1D.m.9     1D.m.7    1D.m.10     1D.m.5     1D.m.6     1D.m.3 
#>  0.0000000  0.2325776  0.4467646  4.8951007 33.5565057 35.1842316 57.4356253 
#>     1D.m.4     1D.m.2     1D.m.1 
#> 57.4356253 57.4356253 57.8942535 
#> 
#> $L.5
#>    2A.m.39    2A.m.40    2A.m.38    2A.m.37    2A.m.36    2A.m.35    2A.m.34 
#>   0.000000   0.000000   0.000000   1.376916   1.836066   1.836066  19.441233 
#>    2A.m.32    2A.m.33    2A.m.30    2A.m.31    2A.m.29    2A.m.28    2A.m.26 
#>  20.108910  20.108910  26.073209  26.073209  41.716710  41.716710  42.634244 
#>    2A.m.27    2A.m.25    2A.m.24    2A.m.23    2A.m.22    2A.m.21    2A.m.20 
#>  42.634244  42.634244  45.848455  49.524783  49.524783  49.524783  54.124872 
#>    2A.m.19    2A.m.18    2A.m.17    2A.m.16    2A.m.15    2A.m.14    2A.m.13 
#>  54.124872  69.248733  69.707163  69.707163  69.707163  69.707163  70.165891 
#>    2A.m.12    2A.m.11    2A.m.10     2A.m.9     2A.m.7     2A.m.8     2A.m.6 
#>  70.624619  71.083348  71.083348  92.475166  92.475166  92.475166  92.961913 
#>     2A.m.4     2A.m.5     2A.m.3     2A.m.2     2A.m.1 
#> 114.318005 114.318005 114.318005 116.545148 123.771568 
#> 
#> $L.6
#>    2B.m.35    2B.m.34    2B.m.33    2B.m.31    2B.m.32    2B.m.30    2B.m.29 
#>  0.0000000  0.4695553  0.4695553  1.8460497  1.8460497  2.3047782  2.7635067 
#>    2B.m.28    2B.m.27    2B.m.26    2B.m.25    2B.m.23    2B.m.24    2B.m.22 
#>  2.7635067  4.1339848  4.5927133  4.5927133  6.5117357  6.6491104  8.3996927 
#>    2B.m.20    2B.m.21    2B.m.19    2B.m.18    2B.m.17    2B.m.16    2B.m.15 
#>  9.8735883 10.0284758 11.7640394 12.6815736 14.0546235 14.0546235 14.5133520 
#>    2B.m.14    2B.m.13    2B.m.12    2B.m.11    2B.m.10     2B.m.9     2B.m.8 
#> 14.5133520 16.3486657 18.1843524 18.1843524 28.4165824 28.4165824 28.4165824 
#>     2B.m.7     2B.m.6     2B.m.5     2B.m.3     2B.m.4     2B.m.1     2B.m.2 
#> 28.8765099 38.1543503 46.9615968 52.0276238 52.0276238 75.6110716 75.7219819 
#> 
#> $L.7
#>   2D1.m.6   2D1.m.7   2D1.m.8   2D1.m.5   2D1.m.4   2D1.m.3   2D1.m.1   2D1.m.2 
#> 0.0000000 0.0000000 0.4587285 4.5664939 5.0280757 5.0280757 5.4903266 5.4903266 
#> 
#> $L.8
#>  2D2.m.13  2D2.m.12  2D2.m.10   2D2.m.9  2D2.m.11   2D2.m.8   2D2.m.7   2D2.m.6 
#>  0.000000  0.000000  6.457687  6.457687  6.457687 36.756312 36.756312 36.756312 
#>   2D2.m.5   2D2.m.4   2D2.m.3   2D2.m.2   2D2.m.1 
#> 38.559460 38.559460 48.313484 48.313484 48.313484 
#> 
#> $L.9
#>    3A.m.37    3A.m.35    3A.m.36    3A.m.33    3A.m.34    3A.m.32    3A.m.31 
#>   0.000000   0.000000   0.000000   1.414973   1.414973   4.131151   4.131151 
#>    3A.m.30    3A.m.29    3A.m.28    3A.m.27    3A.m.26    3A.m.25    3A.m.23 
#>   5.177161   7.525530  25.614865  43.265529  43.265529  49.257350  55.261756 
#>    3A.m.24    3A.m.22    3A.m.21    3A.m.20    3A.m.19    3A.m.16    3A.m.17 
#>  55.261756  62.641921  62.641921  62.641921  63.103611  80.773459  81.232188 
#>    3A.m.18    3A.m.15    3A.m.14    3A.m.13    3A.m.11    3A.m.10     3A.m.9 
#>  81.232188  83.063587  83.063587  83.063587  87.196692  87.196692  87.196692 
#>    3A.m.12     3A.m.7     3A.m.8     3A.m.5     3A.m.6     3A.m.4     3A.m.3 
#>  88.114226 121.204365 121.349262 128.775291 128.973077 133.458782 140.851640 
#>     3A.m.2     3A.m.1 
#> 140.851640 140.851640 
#> 
unlink("MSToutput.txt")
```
