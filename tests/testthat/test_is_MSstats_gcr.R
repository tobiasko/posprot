library(MSstats)

test_that("MSstats demo data set is recognized", {
  QuantData <- suppressMessages(dataProcess(SRMRawData))
  comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0), nrow=1)
  row.names(comparison) <- "T7-T1"
  testResultOneComparison<-suppressMessages(groupComparison(contrast.matrix=comparison, data=QuantData))
  expect_equal(is.MSstats_gcr(testResultOneComparison), TRUE)
  file.remove("msstats.log")
  file.remove("sessionInfo.txt")
})

test_that("my sample data set is recognized",{
  expect_equal(is.MSstats_gcr(sample.data), TRUE)
})

test_that("df is not recognized", {
  df <- data.frame(A = 1:10, B = "Tobi")
  expect_equal(is.MSstats_gcr(df), FALSE)
})
