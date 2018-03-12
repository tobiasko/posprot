library("IRanges")

test_that("helper function finds neoCN pair in IRanges object", {

  ## test instance
  ## ---[R].XXX{E|D}.XXX[R].---

  ## global environment is expected to contain:
  ## test.prot.spec == c("E", "D")
  ## work.prot.spec == c("R")

  test <- IRanges(start = c(1, 6), end = c(5, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label"),
                                         first.aa = c("X", "X"),
                                          last.aa = c("E", "R"),
                                            up.aa = c("R", "E"),
                                          down.aa = c("X", "X")
  )
  expect_is(neoCN(test), "IRanges")
  expect_is(neoCN(test, FALSE), "Pairs")
  expect_equal(start(first(neoCN(test, FALSE))), 6)
  expect_equal(length(neoCN(test)), 1)
  expect_equal(start(neoCN(test)), 1)
  expect_equal(end(neoCN(test)), 10)
})

test_that("helper function finds no pairs", {

  test.prot.spec <- c("G")
  test <- IRanges(start = c(1, 6), end = c(5, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label"),
                                         first.aa = c("X", "X"),
                                          last.aa = c("E", "R"),
                                            up.aa = c("R", "E"),
                                          down.aa = c("X", "X")
  )
  expect_is(neoCN(test), "IRanges")
  #expect_equal(length(.neoCN(test)), 0)
})

test_that("helper function works in absence of known test protease spec", {

  test.prot.spec <- NA
  test <- IRanges(start = c(1, 6), end = c(5, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label"),
                                         first.aa = c("X", "X"),
                                          last.aa = c("E", "R"),
                                            up.aa = c("R", "E"),
                                          down.aa = c("X", "X")
  )

  expect_is(neoCN(test), "IRanges")
  expect_equal(length(neoCN(test)), 1)
  expect_equal(start(neoCN(test)), 1)
  expect_equal(end(neoCN(test)), 10)

})

test_that("main function finds 3 neoCN pairs in list of IRanges", {

  ## list of 3 test instances as above

  test <- IRanges(start = c(1, 6), end = c(5, 10))
  elementMetadata(test) <- data.frame(nterm.label = c("free", "label"),
                                         first.aa = c("X", "X"),
                                          last.aa = c("E", "R"),
                                            up.aa = c("R", "E"),
                                          down.aa = c("X", "X")
  )

  l <- list(a = test, b = test, c = test)
  expect_is(find_neoCN_pairs(l), "list")
  expect_length(find_neoCN_pairs(l), 3)
  expect_is(find_neoCN_pairs(l)[[1]], "IRanges")
  expect_is(find_neoCN_pairs(l, united = FALSE), "list")
  expect_is(find_neoCN_pairs(l, united = FALSE)[[1]], "Pairs")

})
