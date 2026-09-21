test_that("package exports its documented API", {
  exported <- c("mstmap", "mstmap.cross", "mstmap.data.frame", "breakCross",
                "mergeCross", "combineMap", "quickEst", "genClones", "fixClones",
                "heatMap", "pullCross", "pushCross", "subsetCross", "pp.init",
                "pValue", "statGen", "profileGen", "statMark", "profileMark",
                "alignCross")
  for (f in exported) {
    expect_true(exists(f, where = asNamespace("ASMap"), mode = "function"),
                info = paste(f, "is not exported"))
  }
})

test_that("bundled datasets load with the documented shape", {
  expect_s3_class(load_map("mapDH"), "cross")
  expect_s3_class(load_map("mapBC"), "cross")
  expect_s3_class(load_map("mapF2"), "cross")
  expect_s3_class(load_map("mapBCu"), "cross")

  ## mapDHf is the data.frame form: markers in rows, genotypes in columns
  dhf <- load_map("mapDHf")
  expect_s3_class(dhf, "data.frame")
  expect_equal(nrow(dhf), 599)
  expect_equal(ncol(dhf), 218)
})

test_that("mstmap is deterministic", {
  ## The C++ core uses no RNG, so identical input must give identical output.
  ## Every regression test in this suite depends on this property.
  dh <- load_map("mapDH")
  a <- quiet(mstmap.cross(dh, bychr = TRUE, trace = FALSE))
  b <- quiet(mstmap.cross(dh, bychr = TRUE, trace = FALSE))
  expect_identical(map_signature(a), map_signature(b))
})
