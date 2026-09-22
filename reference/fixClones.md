# Consensus genotypes for clonal genotype groups

Consensus genotypes for clonal genotype groups of an R/qtl object.

## Usage

``` r
fixClones(object, gc, id = "Genotype", consensus = TRUE)
```

## Arguments

- object:

  An qtl `object` object with any class structure.

- gc:

  A data frame of genotype clone infomation usually from a call to
  `genClones` (see Details).

- id:

  Character string defining the column of `object$pheno` containing the
  genotype names.

- consensus:

  A logical value. If `TRUE` then consensus genotypes will be calculated
  for each clonal group by intelligently collapsing alleles for each
  marker (see Details). If `FALSE` then for each clonal group the
  genotype with the least missing alleles across the genome will be
  retained and the remaining genotypes from each group will be removed.

## Details

This function provides a very efficient way of dealing with genotype
clones in a genetic marker set. This function can be used at any stage
of the map construction process as it retains linkage group and marker
position information.

The `gc` argument needs to be a data frame of clone information and is
easily obtained from a call to `genClones`. If this function is not used
then the data frame must contain at least three columns with the first
two columns named `"G1"` and `"G2"` containing the pairs of genotypes
that are clones and a `"group"` column that indicates the clonal group
the pairs of genotypes belongs to.

If `consensus = TRUE` then the function will intelligently collapse the
alleles for each marker to form a consensus genotype. Specifically, the
allele value will remain unchanged when there are observed allele values
across all genotypes in the clone group. For cases where there are
missing alleles for some but not all of the genotypes, the consensus
genotype will be given the common allele value from the genotypes that
contained observed allele values. If there is more than one unique
allele value across the genotypes for any marker then it is set to
missing.

## Value

The cross object is returned with identical class structure as the
imputted cross object.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`comparegeno`](https://rdrr.io/pkg/qtl/man/comparegeno.html) and
[`genClones`](https://drj001.github.io/ASMap/reference/genClones.md)

## Examples

``` r

data(mapDH, package = "ASMap")

gc <- genClones(mapDH)
mapDHf <- fixClones(mapDH, gc$cgd, consensus = TRUE)
```
