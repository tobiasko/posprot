test_that("sample data is annotated", {
  n <- get.comparison.labels(sample.data)
  expect_is(annotate_gcr("uniprot.fasta", sample.data, n[1]), "data.frame")
})

test_that("modifications are simplified correctly", {
  expect_equal(simple.mods("ABC[+28.00]DEF"), "ABC[+28]DEF")
  expect_equal(simple.mods("ABC[+28.0]DEF"), "ABC[+28]DEF")
})

test_that("modifications are not simplified", {
  expect_equal(simple.mods("ABC[+28.01]DEF"),  "ABC[+28.01]DEF")
  expect_equal(simple.mods("ABC[+28]DEF"), "ABC[+28]DEF")
})

test_that("modifications are stripped correctly", {
  expect_equal(strip.mods("ABC[+28.01]DEF"), "ABCDEF")
  expect_equal(strip.mods("ABC[-28.01]DEF"), "ABCDEF")
  expect_equal(strip.mods("ABC[-28]DEF"), "ABCDEF")
})

test_that("ends are stripped", {
  expect_equal(strip.ends("_ABCDEF_"), "ABCDEF")
})
