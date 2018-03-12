library(MSstats)

test_that("group comparision results can be extracted from MStstats demo data set", {
  QuantData <- suppressMessages(dataProcess(SRMRawData))
  comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
  row.names(comparison) <- "T7-T1"
  testResultOneComparison<-suppressMessages(groupComparison(contrast.matrix=comparison, data=QuantData))
  expect_is(get.MSstats_gcr(testResultOneComparison, "T7-T1"), "data.frame")
  expect_equal(dim(get.MSstats_gcr(testResultOneComparison, "T7-T1")), c(2, 11))
  file.remove("msstats.log")
  file.remove("sessionInfo.txt")
})
