## Helpers shared across the ASMap test suite.
##
## The linkage map construction in ASMap is deterministic: the C++ MSTmap core
## uses no random number generation, so repeated runs on identical input give
## bitwise identical marker orders and genetic distances. That makes it possible
## to pin current behaviour and detect any change introduced by refactoring.

## mstmap() writes a linkage group summary to stdout from C++ (Rprintf), which
## is not suppressed by suppressMessages(). capture.output() does catch it.
quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- force(expr))
  out
}

## A compact, comparable signature of a constructed linkage map.
##
## Distances are rounded before comparison. The marker ORDER is the primary
## algorithmic output and is compared exactly; the distances are derived
## quantities and are compared with tolerance, because the final cM values are
## floating point reductions that may differ in the last bits between
## compilers and platforms.
map_signature <- function(cross, digits = 6) {
  stopifnot(inherits(cross, "cross"))
  mp <- qtl::pull.map(cross)
  list(
    lg     = names(cross$geno),
    nmar   = as.integer(qtl::nmar(cross)),
    nind   = as.integer(qtl::nind(cross)),
    order  = unname(unlist(lapply(mp, names))),
    dist   = round(unname(unlist(lapply(mp, as.numeric))), digits)
  )
}

## Structural invariants that must hold for any valid linkage map, regardless
## of the particular ordering the algorithm produces. These are the assertions
## that stay meaningful even if the algorithm is deliberately changed.
expect_valid_map <- function(cross) {
  testthat::expect_s3_class(cross, "cross")
  mp <- qtl::pull.map(cross)

  ## every linkage group has at least one marker
  testthat::expect_true(all(qtl::nmar(cross) >= 1))

  ## genetic distances are non-decreasing within each linkage group and start at 0
  for (lg in names(mp)) {
    d <- as.numeric(mp[[lg]])
    testthat::expect_false(any(is.na(d)), info = paste("NA distance in", lg))
    testthat::expect_equal(d[1], 0, info = paste("first marker not at 0 in", lg))
    testthat::expect_true(all(diff(d) >= 0), info = paste("distances not monotonic in", lg))
  }

  ## marker names are unique across the whole map and match the genotype columns
  mn <- unlist(lapply(mp, names))
  testthat::expect_equal(anyDuplicated(mn), 0L)
  for (lg in names(mp)) {
    testthat::expect_identical(names(mp[[lg]]), colnames(cross$geno[[lg]]$data))
  }

  invisible(TRUE)
}

## Datasets used across the suite, loaded on demand.
load_map <- function(name) {
  e <- new.env(parent = emptyenv())
  utils::data(list = name, package = "ASMap", envir = e)
  get(name, envir = e)
}
