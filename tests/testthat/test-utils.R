## Tests for the map manipulation and diagnostic functions.
##
## These are fast and check input validation and structural invariants rather
## than exact numerics, so they are safe to run everywhere.

# ---------------------------------------------------------------------------
# pp.init
# ---------------------------------------------------------------------------

test_that("pp.init returns defaults and validates its arguments", {
  p <- pp.init()
  expect_named(p, c("seg.thresh", "seg.ratio", "miss.thresh", "max.rf", "min.lod"))
  expect_equal(p$seg.thresh, 0.05)
  expect_equal(p$miss.thresh, 0.1)

  expect_error(pp.init(seg.thresh = 1.5), "cannot exceed 1")
  expect_error(pp.init(seg.thresh = "nonsense"), "bonf")
  expect_error(pp.init(miss.thresh = "a"), "wrong type")
  expect_error(pp.init(max.rf = 0.9), "cannot exceeed 0.5")
  expect_error(pp.init(min.lod = 0), "less than or equal to zero")

  ## specifying a ratio clears the threshold
  expect_null(pp.init(seg.ratio = "1:1")$seg.thresh)
  expect_error(pp.init(seg.ratio = "1-1"), "split by")
})

# ---------------------------------------------------------------------------
# breakCross / mergeCross
# ---------------------------------------------------------------------------

test_that("breakCross splits a linkage group at the named marker", {
  dh <- load_map("mapDH")
  lg <- names(nmar(dh))[1]
  mk <- names(pull.map(dh)[[lg]])
  split.at <- mk[10]

  res <- breakCross(dh, split = setNames(list(split.at), lg), suffix = "numeric")

  ## the original group is gone, replaced by two new ones
  expect_false(lg %in% names(nmar(res)))
  expect_true(all(paste(lg, 1:2, sep = ".") %in% names(nmar(res))))

  ## no markers gained or lost
  expect_equal(totmar(res), totmar(dh))
  expect_setequal(markernames(res), markernames(dh))

  ## the split lands in the right place and each new group restarts at 0
  expect_equal(unname(nmar(res)[paste0(lg, ".1")]), 10L)
  for (g in paste(lg, 1:2, sep = ".")) {
    expect_equal(unname(pull.map(res)[[g]][1]), 0)
  }
})

test_that("breakCross validates its input", {
  dh <- load_map("mapDH")
  expect_error(breakCross(dh), "Split cannot be null")
  expect_error(breakCross(dh, split = list(NOSUCH = "m1")), "do not exist")
  expect_error(
    breakCross(dh, split = setNames(list("not-a-marker"), names(nmar(dh))[1])),
    "marker names do not exist")
})

test_that("mergeCross combines linkage groups", {
  dh <- load_map("mapDH")
  g  <- names(nmar(dh))[1:2]
  n_before <- sum(nmar(dh)[g])

  res <- mergeCross(dh, merge = setNames(list(g), "merged"))

  expect_true("merged" %in% names(nmar(res)))
  expect_false(any(g %in% names(nmar(res))))
  expect_equal(unname(nmar(res)["merged"]), n_before)
  expect_equal(totmar(res), totmar(dh))
})

test_that("mergeCross validates its input", {
  dh <- load_map("mapDH")
  expect_error(mergeCross(dh), "Need list")
  expect_error(mergeCross(dh, merge = list(c("1A"))), "2 or greater")
  expect_error(mergeCross(dh, merge = list(c("1A", "NOSUCH"))), "do not appear")
})

test_that("breakCross then mergeCross is marker preserving", {
  dh <- load_map("mapDH")
  lg <- names(nmar(dh))[1]
  mk <- names(pull.map(dh)[[lg]])

  broken <- breakCross(dh, split = setNames(list(mk[10]), lg))
  rejoined <- mergeCross(broken, merge = setNames(list(paste(lg, 1:2, sep = ".")), lg))

  expect_equal(totmar(rejoined), totmar(dh))
  expect_setequal(markernames(rejoined), markernames(dh))
  expect_equal(unname(nmar(rejoined)[lg]), unname(nmar(dh)[lg]))
})

# ---------------------------------------------------------------------------
# subsetCross
# ---------------------------------------------------------------------------

test_that("subsetCross subsets linkage groups and individuals", {
  dh <- load_map("mapDH")
  g  <- names(nmar(dh))[1:3]

  bychr <- subsetCross(dh, chr = g)
  expect_setequal(names(nmar(bychr)), g)
  expect_equal(nind(bychr), nind(dh))

  byind <- subsetCross(dh, ind = 1:50)
  expect_equal(nind(byind), 50)
  expect_equal(totmar(byind), totmar(dh))
})

test_that("subsetCross rejects an invalid ind argument", {
  dh <- load_map("mapDH")
  expect_error(subsetCross(dh, ind = "first"), "logical or numeric")
})

# ---------------------------------------------------------------------------
# statGen / statMark
# ---------------------------------------------------------------------------

test_that("statGen returns per-genotype counts", {
  dh <- load_map("mapDH")
  s  <- statGen(dh, bychr = FALSE, stat.type = c("xo", "dxo", "miss"))

  expect_named(s, c("xo", "dxo", "miss"))
  for (nm in names(s)) {
    expect_length(s[[nm]], nind(dh))
    expect_true(all(s[[nm]] >= 0))
  }
  ## missing counts cannot exceed the number of markers
  expect_true(all(s$miss <= totmar(dh)))
})

test_that("statGen bychr = TRUE returns a matrix per linkage group", {
  dh <- load_map("mapDH")
  s  <- statGen(dh, bychr = TRUE, stat.type = "miss")
  expect_true(is.matrix(s$miss))
  expect_equal(nrow(s$miss), nind(dh))
  expect_equal(ncol(s$miss), nchr(dh))
})

test_that("statGen validates its input", {
  dh <- load_map("mapDH")
  expect_error(statGen(dh, stat.type = "nonsense"), "does not match")
  expect_error(statGen(dh, id = "NoSuchColumn"), "cannot be found")
})

test_that("statMark returns marker and interval statistics", {
  dh <- load_map("mapDH")
  s  <- statMark(dh, stat.type = c("marker", "interval"))

  expect_named(s, c("marker", "interval"))
  expect_equal(nrow(s$marker), totmar(dh))
  expect_true(all(c("chr", "pos", "missing", "dxo") %in% names(s$marker)))
  expect_true(all(s$marker$dxo >= 0))

  expect_true(all(c("erf", "lod", "dist", "mrf", "recomb") %in% names(s$interval)))
  expect_true(all(s$interval$erf >= 0 & s$interval$erf <= 0.5, na.rm = TRUE))
})

test_that("statMark validates its input", {
  dh <- load_map("mapDH")
  expect_error(statMark(dh, stat.type = "nonsense"), "does not match")
})

# ---------------------------------------------------------------------------
# pullCross / pushCross round trip
# ---------------------------------------------------------------------------

test_that("pullCross removes markers and records them", {
  dh <- load_map("mapDH")
  res <- pullCross(dh, type = "co.located")

  expect_true(!is.null(res$co.located))
  expect_true(all(c("table", "data") %in% names(res$co.located)))
  ## pulled markers have left the map
  expect_lt(totmar(res), totmar(dh))
  ## and are accounted for in the stored data
  expect_equal(totmar(res) + ncol(res$co.located$data), totmar(dh))
})

test_that("pullCross then pushCross restores the marker complement", {
  dh <- load_map("mapDH")
  pulled <- pullCross(dh, type = "co.located")
  pushed <- pushCross(pulled, type = "co.located")

  expect_equal(totmar(pushed), totmar(dh))
  expect_setequal(markernames(pushed), markernames(dh))
  ## the co.located slot has been consumed
  expect_null(pushed$co.located)
})

test_that("pullCross with a missing-value threshold selects markers", {
  dh <- load_map("mapDH")
  res <- pullCross(dh, type = "missing", pars = pp.init(miss.thresh = 0.03))
  if (!is.null(res$missing)) {
    expect_lt(totmar(res), totmar(dh))
    expect_equal(totmar(res) + ncol(res$missing$data), totmar(dh))
  } else {
    succeed("no markers exceeded the missing threshold")
  }
})

# ---------------------------------------------------------------------------
# combineMap
# ---------------------------------------------------------------------------

test_that("combineMap merges two maps by genotype", {
  dh <- load_map("mapDH")
  g  <- names(nmar(dh))
  a  <- subsetCross(dh, chr = g[1:3])
  b  <- subsetCross(dh, chr = g[4:6])

  res <- combineMap(a, b, id = "Genotype")

  expect_s3_class(res, "cross")
  expect_equal(totmar(res), totmar(a) + totmar(b))
  expect_setequal(names(nmar(res)), g[1:6])
  expect_equal(nind(res), nind(dh))
})

test_that("combineMap rejects duplicated markers between maps", {
  dh <- load_map("mapDH")
  a  <- subsetCross(dh, chr = names(nmar(dh))[1:2])
  expect_error(combineMap(a, a, id = "Genotype"), "Non-unique markers")
})

# ---------------------------------------------------------------------------
# genClones / fixClones
# ---------------------------------------------------------------------------

test_that("genClones reports genotype similarity", {
  dh <- load_map("mapDH")
  gc <- genClones(dh, tol = 0.95)

  expect_true("cgm" %in% names(gc))
  expect_equal(dim(gc$cgm), c(nind(dh), nind(dh)))

  if (!is.null(gc$cgd)) {
    expect_true(all(c("G1", "G2", "coef", "match", "diff", "group") %in% names(gc$cgd)))
    expect_true(all(gc$cgd$coef > 0.95))
  }
})

test_that("genClones validates its input", {
  dh <- load_map("mapDH")
  expect_error(genClones(dh, id = "NoSuchColumn"), "cannot be found")
  expect_error(genClones(1:10), "class")
})

test_that("fixClones validates its input", {
  dh <- load_map("mapDH")
  expect_error(fixClones(dh), "cannot be missing")
  expect_error(fixClones(dh, gc = 1:10), "data frame")
  expect_error(fixClones(dh, gc = data.frame(a = 1)), "G1")
})

# ---------------------------------------------------------------------------
# quickEst
# ---------------------------------------------------------------------------

test_that("quickEst re-estimates distances without reordering markers", {
  dh <- load_map("mapDH")
  sub <- subsetCross(dh, chr = names(nmar(dh))[1:2])
  est <- quickEst(sub)

  expect_valid_map(est)
  ## marker order is untouched
  expect_identical(markernames(est), markernames(sub))
})

test_that("quickEst rejects non-cross input", {
  expect_error(quickEst(1:10), "class")
})
