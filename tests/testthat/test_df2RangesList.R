test_that("sample data is transformed to RangesList", {
  t <- annotate_gcr("uniprot.fasta", sample.data, label = get.comparison.labels(sample.data))
  expect_is(df2RangesList(t, "uniprot"), "SimpleRangesList")
  expect_equal(IRanges::universe(df2RangesList(t, "uniprot")), "uniprot")
  expect_is(df2RangesList(t)[[1]], "IRanges")
})
