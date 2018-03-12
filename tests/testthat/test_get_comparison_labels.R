library(MSstats)

test_that("MSstats demo data set is recognized", {
  QuantData <- suppressMessages(dataProcess(SRMRawData))
  comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
  row.names(comparison) <- "T7-T1"
  testResultOneComparison<-suppressMessages(groupComparison(contrast.matrix=comparison, data=QuantData))
  expect_equal(get.comparison.labels(testResultOneComparison), "T7-T1")
  file.remove("msstats.log")
  file.remove("sessionInfo.txt")
})

test_that("df generates error",{
  df <- data.frame(A = 1:10, B = "Tobi")
  expect_error(get.comparison.labels(df))
})
