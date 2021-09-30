library(testthat)
library(LDLcalc)


test_that("Test of our Models", {
  stackModel = LDL_ML_Main_StackingAlgorithm(SampleData,0.7,ReportMultiPlot = FALSE,ComparisonPlot = FALSE)
  allModels = LDL_ML_Main_All_Models(SampleData,0.9,ReportMultiPlot = FALSE,ComparisonPlot=FALSE)
  expect_s3_class(stackModel,"caretStack")
  expect_s3_class(allModels,"caretList")

})
