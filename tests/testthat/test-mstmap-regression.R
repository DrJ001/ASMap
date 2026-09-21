## Regression tests for linkage map construction.
##
## These pin the CURRENT output of the MSTmap core so that refactoring work can
## be shown to be behaviour preserving. ASMap is a widely cited package and its
## published results depend on the exact marker orders it produces, so any
## change to those orders must be deliberate and visible.
##
## Snapshot based tests are skip_on_cran(). They compare exact marker orders and
## rounded genetic distances, which are floating point reductions that can
## legitimately differ in the last bits between compilers, optimisation levels
## and platforms. Failing a CRAN check over that would be worse than useless.
## They run locally and in CI, which is where refactoring happens.
##
## The structural invariant tests below have no such sensitivity and run
## everywhere.

# ---------------------------------------------------------------------------
# Exact output pinning (developer regression guard)
# ---------------------------------------------------------------------------

test_that("mstmap.cross output is unchanged for mapDH", {
  skip_on_cran()
  m <- quiet(mstmap.cross(load_map("mapDH"), bychr = TRUE, trace = FALSE))
  expect_snapshot_value(map_signature(m), style = "deparse", tolerance = 1e-6)
})

test_that("mstmap.cross output is unchanged for mapF2 (bcsft/RIL path)", {
  skip_on_cran()
  m <- quiet(mstmap.cross(load_map("mapF2"), bychr = TRUE, trace = FALSE))
  expect_snapshot_value(map_signature(m), style = "deparse", tolerance = 1e-6)
})

test_that("mstmap.cross output is unchanged for mapBC", {
  skip_on_cran()
  m <- quiet(mstmap.cross(load_map("mapBC"), bychr = TRUE, trace = FALSE))
  expect_snapshot_value(map_signature(m), style = "deparse", tolerance = 1e-6)
})

test_that("mstmap.cross output is unchanged for mapBCu (single large group)", {
  skip_on_cran()
  skip_if(Sys.getenv("ASMAP_SLOW_TESTS") == "", "set ASMAP_SLOW_TESTS=1 to run")
  m <- quiet(mstmap.cross(load_map("mapBCu"), bychr = TRUE, trace = FALSE))
  expect_snapshot_value(map_signature(m), style = "deparse", tolerance = 1e-6)
})

test_that("mstmap.data.frame output is unchanged for mapDHf", {
  skip_on_cran()
  m <- quiet(mstmap.data.frame(load_map("mapDHf"), pop.type = "DH",
                               as.cross = TRUE, trace = FALSE))
  expect_snapshot_value(map_signature(m), style = "deparse", tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Structural invariants (platform robust, run everywhere including CRAN)
# ---------------------------------------------------------------------------

test_that("constructed maps satisfy structural invariants", {
  m <- quiet(mstmap.cross(load_map("mapDH"), bychr = TRUE, trace = FALSE))
  expect_valid_map(m)
})

test_that("map construction preserves genotypes and does not invent markers", {
  dh <- load_map("mapDH")
  m  <- quiet(mstmap.cross(dh, bychr = TRUE, trace = FALSE))

  expect_equal(nind(m), nind(dh))

  ## every marker in the output was present in the input
  expect_true(all(markernames(m) %in% markernames(dh)))
  ## with default arguments no markers are dropped
  expect_setequal(markernames(m), markernames(dh))
})

test_that("bychr = FALSE pools markers before clustering", {
  dh <- load_map("mapDH")
  m  <- quiet(mstmap.cross(dh, bychr = FALSE, trace = FALSE))
  expect_valid_map(m)
  ## markers are conserved regardless of how clustering was scoped
  expect_setequal(markernames(m), markernames(dh))
})

# ---------------------------------------------------------------------------
# Relational tests: how outputs must relate to each other under argument
# changes. These express algorithmic intent and do not depend on exact values.
# ---------------------------------------------------------------------------

test_that("haldane and kosambi give the same order but different distances", {
  dh <- load_map("mapDH")
  k <- quiet(mstmap.cross(dh, bychr = TRUE, dist.fun = "kosambi", trace = FALSE))
  h <- quiet(mstmap.cross(dh, bychr = TRUE, dist.fun = "haldane", trace = FALSE))

  ## the distance function rescales the map, it does not reorder it
  expect_identical(map_signature(k)$order, map_signature(h)$order)

  ## but the cM values genuinely differ
  expect_false(isTRUE(all.equal(map_signature(k)$dist, map_signature(h)$dist)))

  ## haldane >= kosambi for the same recombination fraction
  expect_gte(sum(map_signature(h)$dist), sum(map_signature(k)$dist))
})

test_that("objective.fun ML and COUNT both produce valid maps", {
  dh <- load_map("mapDH")
  cnt <- quiet(mstmap.cross(dh, bychr = TRUE, objective.fun = "COUNT", trace = FALSE))
  ml  <- quiet(mstmap.cross(dh, bychr = TRUE, objective.fun = "ML", trace = FALSE))
  expect_valid_map(cnt)
  expect_valid_map(ml)
  expect_setequal(markernames(cnt), markernames(ml))
})

test_that("a stricter p.value splits markers into more linkage groups", {
  dh <- load_map("mapDH")
  loose  <- quiet(mstmap.cross(dh, bychr = FALSE, p.value = 1e-3,  trace = FALSE))
  strict <- quiet(mstmap.cross(dh, bychr = FALSE, p.value = 1e-12, trace = FALSE))
  expect_gte(nchr(strict), nchr(loose))
})

test_that("return.imputed attaches imputed genotypes with consistent dimensions", {
  dh <- load_map("mapDH")
  m  <- quiet(mstmap.cross(dh, bychr = TRUE, return.imputed = TRUE, trace = FALSE))

  expect_true(!is.null(m$imputed.geno))
  expect_setequal(names(m$imputed.geno), names(m$geno))
  for (lg in names(m$imputed.geno)) {
    imp <- m$imputed.geno[[lg]]
    expect_equal(nrow(imp$data), nind(m))
    expect_equal(ncol(imp$data), length(imp$map))
    expect_identical(colnames(imp$data), names(imp$map))
  }
})

test_that("mvest.bc and detectBadData produce valid maps", {
  dh <- load_map("mapDH")
  a <- quiet(mstmap.cross(dh, bychr = TRUE, mvest.bc = TRUE, trace = FALSE))
  b <- quiet(mstmap.cross(dh, bychr = TRUE, detectBadData = TRUE, trace = FALSE))
  expect_valid_map(a)
  expect_valid_map(b)
})

test_that("mstmap.data.frame and mstmap.cross agree on the same data", {
  ## mapDHf is the data.frame form of mapDH
  fromdf <- quiet(mstmap.data.frame(load_map("mapDHf"), pop.type = "DH",
                                    as.cross = TRUE, trace = FALSE))
  expect_valid_map(fromdf)
  expect_setequal(markernames(fromdf), markernames(load_map("mapDH")))
})

test_that("as.cross = FALSE returns a data frame representation", {
  res <- quiet(mstmap.data.frame(load_map("mapDHf"), pop.type = "DH",
                                 as.cross = FALSE, trace = FALSE))
  expect_false(inherits(res, "cross"))
  expect_true(is.data.frame(res$geno))
  expect_true(all(c("markers", "chr", "dist") %in% names(res$geno)))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("mstmap.cross rejects invalid arguments", {
  dh <- load_map("mapDH")
  expect_error(mstmap.cross(dh, dist.fun = "morgan"), "haldane")
  expect_error(mstmap.cross(dh, objective.fun = "MLE"), "COUNT")
  expect_error(mstmap.cross(dh, miss.thresh = 1.5), "between 0 and 1")
  expect_error(mstmap.cross(dh, id = "NoSuchColumn"), "cannot be found")
})

test_that("mstmap.default rejects unsupported classes", {
  expect_error(mstmap(1:10), "data.frame")
})

test_that("mstmap.data.frame rejects invalid arguments", {
  dhf <- load_map("mapDHf")
  expect_error(mstmap.data.frame(dhf, pop.type = "XYZ"), "pop.type|BC")
  expect_error(mstmap.data.frame(dhf, dist.fun = "morgan"), "haldane")
  expect_error(mstmap.data.frame(dhf, objective.fun = "MLE"), "COUNT")
})
