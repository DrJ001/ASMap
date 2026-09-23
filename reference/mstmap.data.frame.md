# Extremely fast linkage map construction for data frame objects using MSTmap.

Extremely fast linkage map construction for data frame objects utilizing
the source code for MSTmap (see Wu et al., 2008). The construction
includes linkage group clustering, marker ordering and genetic distance
calculations.

## Usage

``` r
# S3 method for class 'data.frame'
mstmap(object, pop.type = "DH", dist.fun = "kosambi",
      objective.fun = "COUNT", p.value = 1e-06, noMap.dist = 15,
      noMap.size = 0, miss.thresh = 1, mvest.bc = FALSE,
      detectBadData = FALSE, as.cross = TRUE, return.imputed = FALSE,
      trace = FALSE, ...)
```

## Arguments

- object:

  A `"data.frame"` object containing marker information. The data.frame
  must explicitly be arranged with markers in rows and genotypes in
  columns. Marker names are obtained from the `rownames` of the `object`
  and genotype names are obtained from the `names` component of the
  `object` (see Details).

- pop.type:

  Character string specifying the population type of the data frame
  `object`. Accepted values are `"DH"` (doubled haploid), `"BC"`
  (backcross), `"RILn"` (non-advanced RIL population with n generations
  of selfing) and `"ARIL"` (advanced RIL) (see Details). Default is
  `"DH"`.

- dist.fun:

  Character string defining the distance function used for calculation
  of genetic distances. Options are `"kosambi"` and `"haldane"`. Default
  is `"kosambi"`.

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
  Isolated markers will appear in their own linkage groups ad will be of
  size specified by `noMap.size`.

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
  `"BC","DH","ARIL"` populations only (see Details).

- detectBadData:

  Logical value. If `TRUE` possible genotyping errors are detected, set
  to missing and then imputed as part of the marker ordering algorithm.
  Genotyping errors will also be printed in the file specified by
  `trace`. This is restricted to `"BC","DH","ARIL"` populations only.
  (see Details). Default is `FALSE`.

- as.cross:

  Logical value. If `TRUE` the constructed linkage map is returned as a
  qtl cross object (see Details). If `FALSE` then the constructed
  linkage map is returned as a `data.frame` with extra columns
  indicating the linkage group, marker name/position and genetic
  distance. Default is `TRUE`.

- return.imputed:

  Logical value. If `TRUE` then the imputed marker probability matrix is
  returned for the linkage groups that are constructed (see Details).
  Default is `FALSE`.

- trace:

  An automatic tracing facility. If `trace = FALSE` then minimal
  `MSTmap` output is piped to the screen during the algorithm. If
  `trace = TRUE`, then detailed output from MSTmap is piped to
  "`MSToutput.txt`". This file is equivalent to the output that would be
  obtained from running the MSTmap executable from the command line.

- ...:

  Currently ignored.

## Details

The data frame `object` must have an explicit format with markers in
rows and genotypes in columns. The marker names are required to be in
the `rownames` component and the genotype names are required to be in
the `names` component of the `object`. In each set of names there must
be no spaces. If spaces are detected they are exchanged for a "-". Each
of the columns of the data frame must be of class `"character"` (not
factors). If converting from a matrix, this can easily be achieved by
using the `stringAsFactors = FALSE` argument for any `data.frame`
method.

It is important to know what population type the data frame `object` is
and to correctly input this into `pop.type`. If `pop.type = "ARIL"` then
it is assumed that the minimal number of heterozygotes have been set to
missing before proceeding. The advanced RIL population is then treated
like a backcross population for the purpose of linkage map construction.
Genetic distances are adjusted post construction. For non-advanced RIL
populations `pop.type = "RILn"`, the number of generations of selfing is
limited to 20 to ensure sensible input.

The content of the markers in `object` can either be all numeric (see
below) or all character. If markers are of type character then the
following allelic content must be explicitly adhered to. For `pop.type`
`"BC"`, `"DH"` or `"ARIL"` the two allele types should be represented as
(`"A"` or `"a"`) and (`"B"` or `"b"`). For non-advanced RIL populations
(`pop.type = "RILn"`) phase unknown heterozygotes should be represented
as `"X"`. For all populations, missing marker scores should be
represented as (`"U"` or `"-"`).

This function also extends the functionality of the MSTmap algorithm by
allowing users to input a complete numeric data frame of marker
probabilities for `pop.type` `"BC"`, `"DH"` or `"ARIL"`. The values must
be inclusively between 1 (A) and 0 (B) and be representative of the
probability that the A allele is present. No missing values are allowed.

The algorithm allows an adjustment of the `p.value` threshold for
clustering of markers to distinct linkage groups (see Wu et al., 2008)
and is highly dependent on the number of individuals in the population.
As the number of individuals increases the `p.value` threshold should be
decreased accordingly. This may require some trial and error to achieve
desired results.

If `mvest.bc = TRUE` and the population type is `"BC","DH","ARIL"` then
missing values are imputed before markers are clustered into linkage
groups. This is only a simple imputation that places a 0.5 probability
of the missing observation being one allele or the other and is used to
assist the clustering algorithm when there is known to be high numbers
of missing observations between pairs of markers.

It should be highlighted that for population types `"BC","DH","ARIL"`,
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

For populations `"BC","DH","ARIL"`, if `detectBadData = TRUE`, the
marker ordering algorithm also includes the detection of genotyping
errors. For any individual genotype, the detection method is based on a
weighted Euclidean metric (see Wu et al., 2008) that is a function of
the recombination probabilities of all the markers with the marker
containing the suspicious observation. Any genotyping errors detected
are set to missing and the missing values are then imputed if
`mv.est = TRUE`. Note, the detection of these errors and their amendment
is returned in the imputed probability matrix if `return.imputed = TRUE`

If `as.cross = TRUE` then the constructed object is returned as a qtl
cross object with the appropriate class structure. For `"RILn"`
populations the constructed object is given the class `"bcsft"` by using
the qtl package conversion function `convert2bcsft` with arguments
`F.gen = n` and `BC.gen = 0`. For `"ARIL"` populations the constructed
object is given the class `"riself"`.

If `return.imputed = TRUE` and `pop.type` is one of `"BC","DH","ARIL"`,
then the marker probability matrix is returned for the linkage groups
that have been constructed using the algorithm. Each linkage group is
named identically to the linkage groups of the map and, if
`as.cross = TRUE`, contains an ordered `"map"` element and a `"data"`
element consisting of marker probabilities of the A allele being present
(i.e. P(A) = 1, P(B) = 0). Both elements contain a possibly reduced
version of the marker set that includes all non-colocating markers as
well as the first marker of any set of co-locating markers. If
`as.cross = FALSE` then an ordered data frame of matrix probabilities is
returned.

## Value

If `as.cross = TRUE` the function returns an R/qtl cross object with the
appropriate class structure. The object is a list with usual components
`"pheno"` and `"geno"`. If `as.cross = FALSE` the function returns an
ordered data frame object with additional columns that indicate the
linkage group, the position and marker names and genetic distance of the
markers within in each linkage group. If markers were omitted for any
reason during the construction, the object will have an `"omit"`
component with all omitted markers in a collated matrix. If
`return.imputed = TRUE` then the object will also contain an
`"imputed.geno"` element.

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

[`mstmap.cross`](https://drj001.github.io/ASMap/reference/mstmap.cross.md)

## Examples

``` r

data(mapDH, package = "ASMap")

## forming data frame object from R/qtl object

dfg <- t(do.call("cbind", lapply(mapDH$geno, function(el) el$data)))
dimnames(dfg)[[2]] <- as.character(mapDH$pheno[["Genotype"]])
dfg <- dfg[sample(1:nrow(dfg), nrow(dfg), replace = FALSE),]
dfg[dfg == 1] <- "A"
dfg[dfg == 2] <- "B"
dfg[is.na(dfg)] <- "U"
dfg <- cbind.data.frame(dfg, stringsAsFactors = FALSE)

## construct map

testd <- mstmap(dfg, dist.fun = "kosambi", trace = FALSE)
#> Number of linkage groups: 24
#> The size of the linkage groups are: 33   27  32  33  9   37  4   37  35  56  40  13  41  41  8   6   5   21  54  15  30  10  6   6   
#> The number of bins in each linkage group: 32 27  32  33  8   36  3   36  34  55  40  13  41  41  7   5   4   21  53  14  29  9   5   6   
pull.map(testd)
#> $L1
#>      7A.m.1      7A.m.2      7A.m.3      7A.m.4      7A.m.5      7A.m.6 
#>   0.0000000   0.6767384   2.6942733   7.2943644  15.6275539  20.6906605 
#>      7A.m.7      7A.m.8     7A.m.10      7A.m.9     7A.m.11     7A.m.12 
#>  47.2467766  47.7055051  51.3818405  52.2993747  53.6758691  55.0523636 
#>     7A.m.13     7A.m.14     7A.m.15     7A.m.16     7A.m.18     7A.m.17 
#>  56.0159074  60.5696186  61.9461131  62.4048415  63.7813360  64.2400644 
#>     7A.m.19     7A.m.20     7A.m.21     7A.m.22     7A.m.23     7A.m.24 
#>  64.6987929  65.6163271  87.0175934  96.7745776  97.6940799  99.9892686 
#>     7A.m.25     7A.m.26     7A.m.27     7A.m.28     7A.m.29     7A.m.30 
#> 115.6166115 123.4789448 124.3964789 129.0014952 129.4651077 130.3826419 
#>     7A.m.31     7A.m.33     7A.m.32 
#> 133.8079436 134.4929249 134.4929249 
#> 
#> $L2
#>     6A.m.1     6A.m.2     6A.m.3     6A.m.4     6A.m.5     6A.m.6     6A.m.7 
#>  0.0000000  0.4587285  1.3762626  3.6714514  5.5071381  7.3428248 19.0183998 
#>     6A.m.8     6A.m.9    6A.m.10    6A.m.11    6A.m.12    6A.m.13    6A.m.14 
#> 24.0815064 24.5402349 33.3564951 34.7434160 42.6057493 46.7436102 53.6682798 
#>    6A.m.15    6A.m.16    6A.m.17    6A.m.18    6A.m.19    6A.m.20    6A.m.21 
#> 55.5039665 56.4215007 56.8802292 57.7977634 58.2564918 60.0921785 61.4686730 
#>    6A.m.22    6A.m.23    6A.m.24    6A.m.25    6A.m.26    6A.m.27 
#> 69.3310062 94.1151752 95.0327093 95.4914378 98.2465163 98.7052448 
#> 
#> $L3
#>    4B.m.32    4B.m.31    4B.m.30    4B.m.29    4B.m.28    4B.m.27    4B.m.26 
#>   0.000000   8.333189  21.963866  22.881400  24.257894  29.784884  37.647217 
#>    4B.m.25    4B.m.24    4B.m.23    4B.m.22    4B.m.21    4B.m.20    4B.m.19 
#>  41.785078  42.243806  43.161340  50.086010  51.003544  52.380039  52.838767 
#>    4B.m.18    4B.m.17    4B.m.16    4B.m.15    4B.m.14    4B.m.13    4B.m.12 
#>  53.756301  56.051490  56.510218  58.345905  58.804634  59.263362  59.722091 
#>    4B.m.11    4B.m.10     4B.m.9     4B.m.8     4B.m.7     4B.m.6     4B.m.5 
#>  60.180819  60.639547  61.557082  62.015810  64.310999  78.436651  81.191729 
#>     4B.m.4     4B.m.3     4B.m.2     4B.m.1 
#>  83.946808  91.340123  91.799301 109.472813 
#> 
#> $L4
#>    1B2.m.1    1B2.m.2    1B2.m.3    1B2.m.4    1B2.m.5    1B2.m.6    1B2.m.7 
#>  0.0000000  0.9175342  9.2507237  9.7094521 11.0859466 12.4624410 12.9211695 
#>    1B2.m.8    1B2.m.9   1B2.m.10   1B2.m.11   1B2.m.12   1B2.m.14   1B2.m.13 
#> 14.7568562 16.1333506 16.5920791 17.0508075 17.9683417 20.2635304 21.6400249 
#>   1B2.m.15   1B2.m.16   1B2.m.17   1B2.m.18   1B2.m.19   1B2.m.20   1B2.m.21 
#> 22.0987533 22.5574818 23.9339762 24.3927047 28.0690401 29.9047269 31.7404136 
#>   1B2.m.22   1B2.m.23   1B2.m.24   1B2.m.25   1B2.m.26   1B2.m.27   1B2.m.28 
#> 40.0616812 41.4497817 42.3673159 43.7438103 44.2025388 45.6241701 47.0458014 
#>   1B2.m.29   1B2.m.30   1B2.m.31   1B2.m.32   1B2.m.33 
#> 48.8814881 49.3402166 50.2577508 51.6342452 60.9136461 
#> 
#> $L5
#>   4A.m.9   4A.m.7   4A.m.8   4A.m.6   4A.m.5   4A.m.3   4A.m.4   4A.m.2 
#> 0.000000 1.835687 1.835687 2.294415 2.753144 3.670678 4.129406 6.424595 
#>   4A.m.1 
#> 7.801089 
#> 
#> $L6
#>      3A.m.1      3A.m.2      3A.m.3      3A.m.4      3A.m.5      3A.m.6 
#>   0.0000000   0.4587285   0.9174569   7.8421265  11.5184620  12.4359962 
#>      3A.m.7      3A.m.8     3A.m.12     3A.m.10      3A.m.9     3A.m.11 
#>  18.8936828  19.8112172  50.1097419  53.3251764  54.7016708  56.0781652 
#>     3A.m.13     3A.m.14     3A.m.15     3A.m.17     3A.m.18     3A.m.16 
#>  59.2935997  59.7523281  60.2110566  62.0467433  62.5054718  62.9642003 
#>     3A.m.19     3A.m.20     3A.m.21     3A.m.22     3A.m.23     3A.m.24 
#>  77.0898523  80.3052867  81.2228209  82.1403551  88.5980417  89.0567701 
#>     3A.m.25     3A.m.26     3A.m.27     3A.m.28     3A.m.29     3A.m.30 
#>  95.0485914 101.5062780 107.5156310 119.2094487 129.9208144 138.3605331 
#>     3A.m.31     3A.m.32     3A.m.33     3A.m.34     3A.m.35     3A.m.36 
#> 140.3299205 141.2474547 143.5426434 144.0013724 145.3778673 145.3778673 
#>     3A.m.37 
#> 145.8365957 
#> 
#> $L7
#>   7D.m.1   7D.m.3   7D.m.2   7D.m.4 
#>  0.00000 16.13391 16.13391 21.66090 
#> 
#> $L8
#>     7B.m.1     7B.m.2     7B.m.3     7B.m.4     7B.m.5     7B.m.6     7B.m.7 
#>   0.000000   9.279401  15.271222  15.729951  21.244743  23.552007  24.314460 
#>     7B.m.8     7B.m.9    7B.m.10    7B.m.11    7B.m.12    7B.m.13    7B.m.14 
#>  24.928243  25.845777  26.304506  26.763234  27.680769  28.598303  30.893491 
#>    7B.m.15    7B.m.16    7B.m.17    7B.m.18    7B.m.19    7B.m.20    7B.m.22 
#>  31.352220  32.269754  33.187288  33.646017  35.481703  38.236782  39.154316 
#>    7B.m.21    7B.m.23    7B.m.24    7B.m.25    7B.m.26    7B.m.27    7B.m.28 
#>  39.154316  42.369751  42.828479  45.583558  48.798992  49.716526  68.434452 
#>    7B.m.29    7B.m.30    7B.m.31    7B.m.32    7B.m.33    7B.m.34    7B.m.35 
#>  69.351986  72.107065  78.098886  78.557614  79.016343  80.392837  84.069173 
#>    7B.m.36    7B.m.37 
#>  84.527901 108.162200 
#> 
#> $L9
#>    2B.m.1    2B.m.2    2B.m.3    2B.m.4    2B.m.5    2B.m.6    2B.m.7    2B.m.8 
#>  0.000000  2.779284 25.847324 26.306053 31.833042 41.112443 49.445632 51.281319 
#>    2B.m.9   2B.m.10   2B.m.11   2B.m.12   2B.m.13   2B.m.14   2B.m.15   2B.m.16 
#> 51.740047 52.198776 61.004291 61.921825 63.757512 65.134006 65.592735 66.051463 
#>   2B.m.17   2B.m.18   2B.m.19   2B.m.20   2B.m.21   2B.m.22   2B.m.23   2B.m.24 
#> 66.510191 67.427726 68.345260 69.721754 70.180483 71.556977 72.933472 73.392200 
#>   2B.m.26   2B.m.25   2B.m.27   2B.m.28   2B.m.29   2B.m.30   2B.m.31   2B.m.32 
#> 75.687389 76.146117 76.604846 77.522380 77.981108 78.439837 78.898565 78.898565 
#>   2B.m.34   2B.m.33   2B.m.35 
#> 80.275060 80.733788 97.377493 
#> 
#> $L10
#>      5B.m.1      5B.m.2      5B.m.3      5B.m.4      5B.m.5      5B.m.6 
#>   0.0000000   0.6705738   3.1549806   4.1082794   5.0258136   7.3210023 
#>      5B.m.7      5B.m.8      5B.m.9     5B.m.10     5B.m.11     5B.m.12 
#>  17.0759432  17.5346717  22.1347628  27.1978694  29.9529479  30.8704820 
#>     5B.m.13     5B.m.14     5B.m.15     5B.m.16     5B.m.17     5B.m.18 
#>  44.0087763  44.4675048  59.0908447  59.0908447  59.5495731  63.2259086 
#>     5B.m.19     5B.m.20     5B.m.21     5B.m.22     5B.m.23     5B.m.24 
#>  64.1830094  66.0883037  67.4947790  70.7102134  76.2372026  77.6136970 
#>     5B.m.25     5B.m.26     5B.m.27     5B.m.28     5B.m.29     5B.m.30 
#>  78.0724255  78.5311539  82.6690149 105.7370555 106.6545896 107.1133181 
#>     5B.m.31     5B.m.32     5B.m.33     5B.m.34     5B.m.35     5B.m.36 
#> 108.0308523 108.9483864 109.4071149 110.3246491 110.7833775 111.2421060 
#>     5B.m.37     5B.m.38     5B.m.39     5B.m.40     5B.m.41     5B.m.42 
#> 113.5374089 114.0207066 114.9627093 115.4214445 115.8801730 116.3389015 
#>     5B.m.43     5B.m.44     5B.m.45     5B.m.46     5B.m.47     5B.m.48 
#> 116.7976299 117.7151641 118.1738926 118.6326210 119.0913495 119.5500780 
#>     5B.m.49     5B.m.50     5B.m.51     5B.m.52     5B.m.53     5B.m.54 
#> 120.0088064 120.4675349 122.3032216 122.7619501 125.0571388 128.7334743 
#>     5B.m.55     5B.m.56 
#> 141.8717685 142.7893027 
#> 
#> $L11
#>    2A.m.1    2A.m.2    2A.m.3    2A.m.4    2A.m.5    2A.m.6    2A.m.7    2A.m.8 
#>   0.00000  14.62334  24.37828  25.29581  27.13150  41.25715  46.78414  47.24287 
#>    2A.m.9   2A.m.10   2A.m.11   2A.m.12   2A.m.13   2A.m.16   2A.m.15   2A.m.14 
#>  51.38073  67.51465  67.97337  68.43210  68.89083  71.64591  74.40099  74.87506 
#>   2A.m.17   2A.m.18   2A.m.19   2A.m.20   2A.m.21   2A.m.22   2A.m.23   2A.m.24 
#>  75.33379  77.16948  89.81787  90.27660  94.41446  94.87319  95.33192  98.54735 
#>   2A.m.25   2A.m.26   2A.m.27   2A.m.28   2A.m.29   2A.m.30   2A.m.31   2A.m.32 
#> 100.84254 103.59762 105.43331 106.35084 107.29001 122.91735 123.37608 128.43918 
#>   2A.m.33   2A.m.34   2A.m.35   2A.m.36   2A.m.37   2A.m.38   2A.m.39   2A.m.40 
#> 128.89791 134.35018 145.98741 147.36390 148.28254 149.66120 150.58250 153.23654 
#> 
#> $L12
#>    2D2.m.1    2D2.m.2    2D2.m.3    2D2.m.4    2D2.m.5    2D2.m.6    2D2.m.7 
#>  0.0000000  0.4587285  0.9218419  8.7809685 11.0774159 13.3753890 13.8368962 
#>    2D2.m.8    2D2.m.9   2D2.m.10   2D2.m.11   2D2.m.12   2D2.m.13 
#> 17.9787481 49.5879130 50.0466418 50.5081742 57.4299864 60.6454208 
#> 
#> $L13
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
#> $L14
#>     6B.m.1     6B.m.2     6B.m.3     6B.m.4     6B.m.5     6B.m.6     6B.m.7 
#>   0.000000   0.918229   1.377097   3.212784   4.589279   5.048007  14.327408 
#>     6B.m.8     6B.m.9    6B.m.10    6B.m.11    6B.m.12    6B.m.13    6B.m.14 
#>  34.105598  34.564326  35.023055  41.947724  46.085585  47.003119  47.461848 
#>    6B.m.15    6B.m.16    6B.m.20    6B.m.19    6B.m.18    6B.m.17    6B.m.22 
#>  48.838342  49.755876  51.132371  52.618139  53.426456  53.885185  55.261679 
#>    6B.m.21    6B.m.23    6B.m.24    6B.m.25    6B.m.26    6B.m.27    6B.m.28 
#>  55.720408  56.637942  59.393020  59.851749  67.710513  68.172724  68.631453 
#>    6B.m.29    6B.m.30    6B.m.31    6B.m.32    6B.m.33    6B.m.34    6B.m.35 
#>  69.090181  71.845260  89.002117  92.252955  94.583475 101.976331 104.271519 
#>    6B.m.36    6B.m.37    6B.m.38    6B.m.39    6B.m.40    6B.m.41 
#> 104.730248 105.688281 107.943156 110.717806 112.149432 115.411344 
#> 
#> $L15
#>  2D1.m.8  2D1.m.7  2D1.m.6  2D1.m.5  2D1.m.4  2D1.m.3  2D1.m.2  2D1.m.1 
#>  0.00000 13.13829 13.13829 16.76419 17.73191 18.22165 18.72569 19.20095 
#> 
#> $L16
#>     5D.m.1     5D.m.2     5D.m.3     5D.m.5     5D.m.4     5D.m.6 
#>  0.0000000  0.4587285  0.9174569 17.5611621 17.5611621 18.0198906 
#> 
#> $L17
#>  1B1.m.4  1B1.m.5  1B1.m.3  1B1.m.2  1B1.m.1 
#> 0.000000 0.000000 4.137861 6.892939 7.810474 
#> 
#> $L18
#>  4A.m.30  4A.m.29  4A.m.28  4A.m.27  4A.m.26  4A.m.25  4A.m.24  4A.m.23 
#>  0.00000  6.92467 18.60024 26.93343 43.06735 44.44384 44.90257 45.82010 
#>  4A.m.22  4A.m.21  4A.m.20  4A.m.19  4A.m.18  4A.m.17  4A.m.16  4A.m.15 
#> 46.73764 47.19637 48.11390 51.79024 53.16673 53.62546 83.92398 84.38271 
#>  4A.m.14  4A.m.13  4A.m.12  4A.m.11  4A.m.10 
#> 86.67790 94.54023 94.99896 95.91650 98.21169 
#> 
#> $L19
#>     5A.m.1     5A.m.2     5A.m.3     5A.m.4     5A.m.5     5A.m.6     5A.m.7 
#>   0.000000   6.457687   6.916415   8.752102   9.210830   9.669559  10.587093 
#>     5A.m.8     5A.m.9    5A.m.10    5A.m.11    5A.m.12    5A.m.13    5A.m.14 
#>  15.650199  21.642021  51.296225  51.754954  53.590640  58.190732  62.328593 
#>    5A.m.15    5A.m.16    5A.m.17    5A.m.18    5A.m.19    5A.m.20    5A.m.21 
#>  64.164279  72.497469  73.873963  75.250458  75.709186  76.167915  81.231021 
#>    5A.m.22    5A.m.23    5A.m.24    5A.m.25    5A.m.26    5A.m.28    5A.m.27 
#>  94.861697  95.779232 101.771053 102.229781 105.906117 108.201306 108.660034 
#>    5A.m.29    5A.m.30    5A.m.31    5A.m.32    5A.m.33    5A.m.34    5A.m.37 
#> 109.577568 111.413255 111.871983 112.789518 115.544596 118.760030 120.136525 
#>    5A.m.36    5A.m.35    5A.m.38    5A.m.40    5A.m.39    5A.m.41    5A.m.42 
#> 120.136525 120.595253 121.512788 123.807976 124.266705 128.866796 130.243290 
#>    5A.m.43    5A.m.44    5A.m.45    5A.m.46    5A.m.47    5A.m.48    5A.m.49 
#> 130.702019 131.160747 131.619476 132.078204 140.883719 142.719406 143.178134 
#>    5A.m.50    5A.m.51    5A.m.52    5A.m.53    5A.m.54 
#> 162.424129 180.097313 181.017554 181.476283 181.935011 
#> 
#> $L20
#>     6D.m.1     6D.m.2     6D.m.3     6D.m.5     6D.m.4     6D.m.6     6D.m.7 
#>  0.0000000  0.4587285  0.4587285  1.3762626  1.8349911  2.2937196  2.7524480 
#>     6D.m.8     6D.m.9    6D.m.10    6D.m.11    6D.m.12    6D.m.13    6D.m.14 
#> 12.9853614 13.4447453 14.3622795 32.6007449 61.5944625 68.0736298 71.2890642 
#>    6D.m.15 
#> 72.2115715 
#> 
#> $L21
#>   3B.m.30   3B.m.29   3B.m.28   3B.m.27   3B.m.26   3B.m.25   3B.m.24   3B.m.23 
#>   0.00000  10.23223  22.39308  30.25542  34.39367  34.85278  38.99064  53.45212 
#>   3B.m.22   3B.m.21   3B.m.20   3B.m.19   3B.m.17   3B.m.18   3B.m.16   3B.m.15 
#>  56.35658  61.41969  66.48279  67.40033  68.31786  68.77659  69.23532  69.69405 
#>   3B.m.14   3B.m.13   3B.m.12   3B.m.11   3B.m.10    3B.m.9    3B.m.8    3B.m.7 
#>  73.83191  74.29064  74.74937  74.74937  75.20810  97.15983  98.07737 104.62584 
#>    3B.m.6    3B.m.5    3B.m.4    3B.m.3    3B.m.2    3B.m.1 
#> 129.29707 131.59226 132.05099 132.50972 133.42725 133.88598 
#> 
#> $L22
#>  1D.m.10   1D.m.9   1D.m.8   1D.m.7   1D.m.6   1D.m.5   1D.m.3   1D.m.4 
#>  0.00000 22.50673 22.96620 23.88404 55.53218 59.75146 78.46976 78.46976 
#>   1D.m.2   1D.m.1 
#> 78.92881 81.22400 
#> 
#> $L23
#>    3D.m.2    3D.m.1    3D.m.3    3D.m.4    3D.m.5    3D.m.6 
#> 0.0000000 0.0000000 0.4606647 3.6788493 8.2797557 8.7384897 
#> 
#> $L24
#>    4D.m.1    4D.m.2    4D.m.3    4D.m.4    4D.m.5    4D.m.6 
#>  0.000000  1.092830  2.600342  4.480309  8.662696 13.262787 
#> 

## let's get a timing on that ...

system.time(testd <- mstmap(dfg, dist.fun = "kosambi", trace = FALSE))
#> Number of linkage groups: 24
#> The size of the linkage groups are: 33   27  32  33  9   37  4   37  35  56  40  13  41  41  8   6   5   21  54  15  30  10  6   6   
#> The number of bins in each linkage group: 32 27  32  33  8   36  3   36  34  55  40  13  41  41  7   5   4   21  53  14  29  9   5   6   
#>    user  system elapsed 
#>   0.433   0.000   0.433 
```
