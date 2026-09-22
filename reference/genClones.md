# Find and report genotype clones

Find and report genotype clones for qtl objects.

## Usage

``` r
genClones(object, chr, tol = 0.9, id = "Genotype")
```

## Arguments

- object:

  An qtl `object` object with any class structure.

- chr:

  A character string of linkage group names.

- tol:

  Pairs of genotypes with a proporion of matching alleles above this
  tolerance will be returned.

- id:

  Character string defining the column of `object$pheno` containing the
  genotype names.

## Details

This function extends the functionality of `comparegeno` in the qtl
package by providing breakdown statistics for the pairs of genotypes
that have a proportion of matching alleles above `tol`.

## Value

A list is returned with the matrix from `comparegeno` as an element
`cgm` and the breakdown statistics for returned genotype pairs in `cgd`.
Specifically, the statistics contain a `"group"` column which determines
the clonal group the pair of genotypes belongs to.

## References

Taylor, J., Butler, D. (2017) R Package ASMap: Efficient Genetic Linkage
Map Construction and Diagnosis. Journal of Statistical Software,
**79**(6), 1–29.

## Author

Julian Taylor

## See also

[`comparegeno`](https://rdrr.io/pkg/qtl/man/comparegeno.html) and
[`fixClones`](https://drj001.github.io/ASMap/reference/fixClones.md)

## Examples

``` r

data(mapDH, package = "ASMap")

gc <- genClones(mapDH)
```
