## Tests for defects found during code review and fixed in 1.1-0.
##
## Each test asserts the CORRECTED behaviour. They previously lived in
## test-known-defects.R, where they documented the broken behaviour, so that
## fixing each defect produced a visible failing test.

# ---------------------------------------------------------------------------
# `type` resolved before first use (match.arg)
# ---------------------------------------------------------------------------

test_that("pushCross() works with its default type argument", {
  ## Previously `if (type != "unlinked")` was evaluated while type was still
  ## the length-4 default vector. R >= 4.2 raises "the condition has length > 1",
  ## so pushCross() failed outright for anyone relying on the default.
  dh <- load_map("mapDH")
  pulled <- pullCross(dh, type = "co.located")

  expect_no_error(pushed <- pushCross(pulled))

  ## the default must behave exactly as the first choice, "co.located"
  explicit <- pushCross(pulled, type = "co.located")
  expect_equal(totmar(pushed), totmar(explicit))
  expect_setequal(markernames(pushed), markernames(explicit))

  ## and it round trips back to the original marker complement
  expect_setequal(markernames(pushed), markernames(dh))
})

test_that("pullCross() works with its default type argument", {
  ## Previously object[[type]] was evaluated with the length-3 default vector,
  ## which performs recursive indexing and failed with "no such index at level 1".
  dh <- load_map("mapDH")

  expect_no_error(pulled <- pullCross(dh))

  explicit <- pullCross(dh, type = "co.located")
  expect_equal(totmar(pulled), totmar(explicit))
  expect_setequal(markernames(pulled), markernames(explicit))
})

test_that("pushCross() reports clearly when there is nothing to push back", {
  ## The genotype-count check dereferenced object[[type]] before establishing
  ## that it exists. dim(NULL)[1] is NULL, so the comparison produced
  ## logical(0) and "argument is of length zero", masking the intended message.
  dh <- load_map("mapDH")
  expect_error(pushCross(dh), "no markers of this type to push back")
  expect_error(pushCross(dh, type = "seg.distortion"), "no markers of this type to push back")
})

test_that("pushCross() and pullCross() still accept partial matching", {
  dh <- load_map("mapDH")
  pulled <- pullCross(dh, type = "co")
  expect_setequal(markernames(pulled), markernames(pullCross(dh, type = "co.located")))
})

# ---------------------------------------------------------------------------
# Spaces in marker and genotype names
# ---------------------------------------------------------------------------

test_that("spaces in marker names produce a warning, not an error", {
  ## The warning() call contained unescaped quotes, which R parsed as a
  ## subtraction of two strings, failing with "non-numeric argument to binary
  ## operator" whenever the branch was reached.
  dhf <- load_map("mapDHf")
  bad <- dhf
  rownames(bad)[3] <- "M 3 with spaces"

  expect_warning(
    res <- quiet(mstmap.data.frame(bad, pop.type = "DH", as.cross = TRUE, trace = FALSE)),
    "Replacing spaces in marker names")

  ## the space has been replaced and the marker survives into the map
  expect_true("M-3-with-spaces" %in% markernames(res))
  expect_false(any(grepl(" ", markernames(res))))
})

test_that("spaces in genotype names produce a warning, not an error", {
  ## The corrected names were assigned to rownames(object) rather than
  ## names(object). Markers are rows and genotypes are columns, so the lengths
  ## differed and the assignment failed with "invalid 'row.names' length".
  dhf <- load_map("mapDHf")
  bad <- dhf
  names(bad)[5] <- "G 5 with spaces"

  expect_warning(
    res <- quiet(mstmap.data.frame(bad, pop.type = "DH", as.cross = TRUE, trace = FALSE)),
    "Replacing spaces in genotype names")

  ## the genotype name is corrected, and marker names are left alone
  expect_true("G-5-with-spaces" %in% as.character(res$pheno$Genotype))
  expect_setequal(markernames(res), markernames(load_map("mapDH")))
})

test_that("correcting names does not change the resulting map", {
  ## Renaming a genotype must not alter linkage map construction.
  dhf <- load_map("mapDHf")
  bad <- dhf
  names(bad)[5] <- "G 5 with spaces"

  base  <- quiet(mstmap.data.frame(dhf, pop.type = "DH", as.cross = TRUE, trace = FALSE))
  fixed <- suppressWarnings(
    quiet(mstmap.data.frame(bad, pop.type = "DH", as.cross = TRUE, trace = FALSE)))

  expect_identical(map_signature(base)$order, map_signature(fixed)$order)
  expect_equal(map_signature(base)$dist, map_signature(fixed)$dist)
})

# ---------------------------------------------------------------------------
# alignCross argument checking
# ---------------------------------------------------------------------------

test_that("alignCross() reports a helpful error for a non-list maps argument", {
  ## if(!grep(...)) raised "argument is of length zero" when there was no match,
  ## masking the intended message.
  dh <- load_map("mapDH")
  expect_error(alignCross(dh, maps = dh), "must be a list of maps")
})

test_that("alignCross() rejects a non-cross object without an NA condition", {
  ## class(object)[2] != "cross" is NA for objects with a single class,
  ## producing "missing value where TRUE/FALSE needed".
  expect_error(alignCross(1:10, maps = list(a = 1)), "should have class")
})

test_that("alignCross() accepts a list of reference maps", {
  dh <- load_map("mapDH")
  ref <- subsetCross(dh, chr = names(nmar(dh))[1:3])
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(res <- alignCross(dh, maps = list(ref = ref)))
  expect_true(all(c("map", "marker", "ref.chr", "ref.dist") %in% names(res)))
})

# ---------------------------------------------------------------------------
# Miscellaneous
# ---------------------------------------------------------------------------

test_that("pValue() legend labels are correctly separated", {
  ## paste(dist, " cM", sepm = "") passed sepm through ... instead of setting
  ## the separator, giving "25  cM " rather than "25 cM".
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(pValue(dist = c(25, 30), pop.size = 100:150))
  expect_equal(paste(25, " cM", sep = ""), "25 cM")
})

test_that("statMark() error message points at the right help page", {
  dh <- load_map("mapDH")
  expect_error(statMark(dh, stat.type = "nonsense"), "\\?statMark")
})

test_that("mstmap.data.frame() rejects numeric input for RILn populations", {
  ## The guard used is.numeric() on a data.frame, which is always FALSE, so the
  ## intended error could never be raised.
  num <- as.data.frame(matrix(runif(20), nrow = 5))
  expect_error(mstmap.data.frame(num, pop.type = "RIL2"),
               "not available for RILn")
})

test_that("mstmap methods do not pick up a stray global omit.list", {
  ## exists("omit.list") searched enclosing environments, so an unrelated
  ## object of that name in the global environment would be used.
  omit.list <- "stray global"
  on.exit(rm(omit.list), add = TRUE)

  dh <- load_map("mapDH")
  res <- quiet(mstmap.cross(dh, bychr = TRUE, trace = FALSE))
  expect_false(identical(res$omit, "stray global"))
})
