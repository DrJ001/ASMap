## Known defects.
##
## These tests document CURRENT, INCORRECT behaviour. They are here so that
## fixing each defect produces a visible, failing test that must be updated to
## assert the corrected behaviour. Each block records what the code does today
## and, in a comment, what it should do instead.
##
## When a defect is fixed, replace the expectation with the "should be" form
## and move the test into the appropriate file.

test_that("DEFECT: pushCross() uses `type` before match.arg()", {
  ## `if (type != "unlinked")` is evaluated while `type` is still the full
  ## default vector of length 4. R >= 4.2 makes a length > 1 condition an
  ## error, so calling pushCross() with the default argument fails outright.
  ##
  ## SHOULD BE: type <- match.arg(type) before first use, so that
  ## pushCross(pulled) behaves as pushCross(pulled, type = "co.located").
  dh <- load_map("mapDH")
  pulled <- pullCross(dh, type = "co.located")
  expect_error(pushCross(pulled), "the condition has length > 1")
})

test_that("DEFECT: pullCross() indexes with the whole default `type` vector", {
  ## object[[type]] with type of length 3 performs recursive indexing rather
  ## than selecting a single component, because match.arg() has not run yet.
  ##
  ## SHOULD BE: type <- match.arg(type) before object[[type]] is touched.
  dh <- load_map("mapDH")
  expect_error(pullCross(dh), "no such index|subscript out of bounds")
})

test_that("DEFECT: alignCross() uses if(!grep(...)) on a possibly empty match", {
  ## grep() returns integer(0) when there is no match, and if(!integer(0))
  ## raises "argument is of length zero" instead of the intended message.
  ##
  ## SHOULD BE: if (!length(grep("list", mp))) stop("maps argument must be a list of maps")
  dh <- load_map("mapDH")
  expect_error(alignCross(dh, maps = dh), "argument is of length zero")
})

test_that("DEFECT: mstmap.data.frame() errors on spaces in marker names", {
  ## The intended behaviour is to replace spaces with "-" and warn. The
  ## warning() call contains unescaped quotes, so R parses the message as
  ## a subtraction of two strings, which fails at run time.
  ##
  ## SHOULD BE: a warning, and construction continues with corrected names.
  dhf <- load_map("mapDHf")
  bad <- dhf
  rownames(bad)[3] <- paste("M 3 with spaces")
  expect_error(
    suppressWarnings(mstmap.data.frame(bad, pop.type = "DH", trace = FALSE)),
    "non-numeric argument to binary operator")
})

test_that("DEFECT: mstmap.data.frame() assigns genotype names to rownames", {
  ## The block handling spaces in genotype names writes the corrected names to
  ## rownames(object) rather than names(object). Because the data frame has
  ## markers in rows and genotypes in columns, the lengths differ and the
  ## assignment fails.
  ##
  ## SHOULD BE: names(object) <- gsub(" ", "-", names(object))
  dhf <- load_map("mapDHf")
  bad <- dhf
  names(bad)[5] <- "G 5 with spaces"
  expect_error(
    suppressWarnings(mstmap.data.frame(bad, pop.type = "DH", trace = FALSE)),
    "invalid 'row.names' length")
})

test_that("DEFECT: pValue() legend labels contain a paste() typo", {
  ## paste(dist, " cM", sepm = "") passes `sepm` through ... instead of
  ## setting the separator, so labels gain stray spaces.
  ##
  ## SHOULD BE: paste(dist, " cM", sep = "") giving "25 cM".
  expect_equal(paste(25, " cM", sepm = ""), "25  cM ")
  expect_equal(paste(25, " cM", sep = ""), "25 cM")
})

test_that("DEFECT: is.numeric() on a data.frame makes the RILn guard dead code", {
  ## mstmap.data.frame() guards the RILn branch with if (is.numeric(object)),
  ## but object is always a data.frame, for which is.numeric() is always FALSE.
  ## The intended error can therefore never be raised.
  ##
  ## SHOULD BE: if (all(sapply(object, is.numeric)))
  df <- data.frame(a = 1:3, b = 4:6)
  expect_false(is.numeric(df))
  expect_true(all(sapply(df, is.numeric)))
})
