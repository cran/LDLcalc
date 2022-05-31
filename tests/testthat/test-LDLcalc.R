library(testthat)
library(LDLcalc)
library(data.table)

test_that("Test All9Models & Stacking Algorithm Model", {
  stackModel = LDL_ML_Main_StackingAlgorithm(SampleData,0.7,ReportMultiPlot = FALSE,ComparisonPlot = FALSE)
  allModels = LDL_ML_Main_All_Models(SampleData,0.9,ReportMultiPlot = FALSE,ComparisonPlot=FALSE)
  expect_s3_class(stackModel,"caretStack")
  expect_s3_class(allModels,"caretList")

})

test_that("Test 9ML_models_functions", {
  modelPrediction = LDL_ML_predict(model$model,data.table::data.table(CHOL=170.5,HDL=35.12,TG=175))
  expect_length(modelPrediction,1)

  model = LDL_ML_Main(SampleData,0.7,"lm",ReportMultiPlot=FALSE)
  expect_type(model,"list")

})


test_that("Test Aip.R function", {
  AIPcalculation = AIPcalc(sampleA$TG,sampleA$HDL)
  expect_length(AIPcalculation,34)

  AIPpropagationVar = AIPErrPrp(sampleA$TG,sampleA$HDL)
  expect_length(AIPpropagationVar,1)

  AIPpropagationVarTaylor = AIPErrPrp2Ord(sampleA$TG,sampleA$HDL)
  expect_length(AIPpropagationVar,1)

  AIPbootstrVar = AIPbootVrnc(sampleA$TG,sampleA$HDL,noOfReps = 2)
  expect_type(AIPbootstrVar,"list")


})



test_that("Test of AIPerrorVariance.R function", {
  dfHDL = CV_Range(sampleB$HDL,0,30,maxRandIter = 10, plot=FALSE)
  expect_s3_class(dfHDL,"data.frame")

  AIP_HDLVrnc=AIP_HDLVrnc(dfHDL,sampleA$TG,bootStrpReps=2)
  expect_type(AIP_HDLVrnc,"list")

  dfTG = CV_Range(sampleA$TG,0,30,maxRandIter = 2, plot=FALSE)
  expect_s3_class(dfTG,"data.frame")

  AIP_TGVrnc=AIP_TGVrnc(dfTG,sampleA$HDL,bootStrpReps=2)
  expect_type(AIP_TGVrnc,"list")
})


test_that("Test CVfunctionRange function", {
  DataFrame=CV_Range(sampleA$LDL,0,10,maxRandIter = 2, plot=FALSE)
  expect_s3_class(DataFrame,"data.frame")


  CV=CV(sampleA$LDL)
  expect_length(CV,1)
})




test_that("Test fErrorCorCov function", {
  ErrorOFCorCov=ErrorCorCov(sampleA$HDL[1:20],sampleA$CHOL[1:20],plot = FALSE)
  expect_type(ErrorOFCorCov,"list")


})


test_that("Test Jensen-Shannon.R function", {
  JSDCalc = JSD(model[["trainData"]]$LDLd,model[["testData"]]$LDLd)
  expect_type(JSDCalc,"list")

  JSD.between.empirical.Normal = JSDNormal(sampleA,"LDL")
  expect_type(JSD.between.empirical.Normal,"list")

  chebysevBounds = chebyshev(sampleA$LDL)
  expect_type(chebysevBounds,"list")



})


test_that("Test LDLErrorPropagation&boostrapVariance function", {
  LDLboostrpVar = LDLbootVrnc(sampleA$CHOL,sampleA$HDL,sampleA$TG,noOfReps = 2)
  expect_type(LDLboostrpVar,"list")

  LDLerrorPrp = LDLErrPrp(sampleA$CHOL,sampleA$HDL,sampleA$TG)
  expect_length(LDLerrorPrp,1)

  dfCHOL = CV_Range(sampleA$CHOL,0,30,maxRandIter = 2, plot=FALSE)
  expect_s3_class(dfCHOL,"data.frame")

  LDLCHOLVar = LDL_CHOLVrnc(dfCHOL,sampleA$HDL,sampleA$TG,bootStrpReps=2)
  expect_type(LDLCHOLVar,"list")

  dfTG = CV_Range(sampleA$TG,0,30,maxRandIter = 2, plot=FALSE)
  expect_s3_class(dfTG,"data.frame")

  LDLTGVar=LDL_TGVrnc(dfTG,sampleA$CHOL,sampleA$HDL,bootStrpReps=2)
  expect_type(LDLTGVar,"list")

  dfHDL = CV_Range(sampleA$HDL,0,30,maxRandIter = 10, plot=FALSE)
  expect_s3_class(dfHDL,"data.frame")

  LDLHDLVar=LDL_HDLVrnc(dfHDL,sampleA$CHOL,sampleA$TG,bootStrpReps=2)
  expect_type(LDLHDLVar,"list")

})

test_that("plotK.R", {
  TGVariances = CV_Range(sampleA$TG,15,16,plot=FALSE)
  LDLTGSampleDependance = LDL_TGVrnc(TGVariances,sampleA$CHOL, sampleA$HDL, bootStrpReps =2)
  expect_type(LDLTGSampleDependance,"list")
  PLTTGvrncLDL=plotTGVrncToLDL(TGVariances,LDLTGSampleDependance$ErrPropVrnc,LDLTGSampleDependance$BootVrnc)
  expect_type(PLTTGvrncLDL,"list")

  TGVariances = CV_Range(sampleA$TG,15,16,plot=FALSE)
  AIPVrncChngTGVrnc = AIP_TGVrnc(TGVariances,sampleA$HDL,bootStrpReps = 2)
  TGErrPropVrnc = AIPVrncChngTGVrnc$ErrPropVrnc
  TGErrPropVrnc2Ord = AIPVrncChngTGVrnc$ErrPropVrnc2Ord
  TGBootVrnc = AIPVrncChngTGVrnc$BootVrnc
  PLTAIPvrncTG=plotAIP_TGVrnc(TGVariances,TGErrPropVrnc,TGErrPropVrnc2Ord,TGBootVrnc)
  expect_type(PLTAIPvrncTG,"list")


  CHOLVariances = CV_Range(sampleA$CHOL,1,10,plot=FALSE)
  LDLCHOLDependance = LDL_CHOLVrnc(CHOLVariances, sampleA$HDL, sampleA$TG, bootStrpReps=2)
  PLTCHOLvrncLDL=plotCholVrncToLDL(CHOLVariances,LDLCHOLDependance$ErrPropVrnc,LDLCHOLDependance$BootVrnc)
  expect_type(PLTCHOLvrncLDL,"list")

  HDLVariances = CV_Range(sampleA$HDL,15,25,plot=FALSE)
  LDLHDLSampleDependance = LDL_HDLVrnc(HDLVariances,sampleA$CHOL, sampleA$TG, bootStrpReps=2)
  PLTHDLvrncLDL=plotHDLVrncToLDL(HDLVariances,LDLHDLSampleDependance$ErrPropVrnc,LDLHDLSampleDependance$BootVrnc)
  expect_type(PLTHDLvrncLDL,"list")


  AIPVrncChngHDLVrnc = AIP_HDLVrnc(HDLVariances,sampleA$TG, bootStrpReps=2)
  HDLErrPropVrnc = AIPVrncChngHDLVrnc$ErrPropVrnc
  HDLErrPropVrnc2Ord = AIPVrncChngHDLVrnc$ErrPropVrnc2Ord
  HDLBootVrnc = AIPVrncChngHDLVrnc$BootVrnc
  PLTAIPvrcHDL=plotAIP_HDLVrnc(HDLVariances,HDLErrPropVrnc,HDLErrPropVrnc2Ord,HDLBootVrnc)
  expect_type(PLTAIPvrcHDL,"list")

})

test_that("Plots.R", {
  LDLbootstrVar=as.data.frame(LDLbootVrnc(sampleA$CHOL,sampleA$HDL,sampleA$TG,noOfReps = 2))
  PLTdensOfVar=DensityPlotOfVar(LDLbootstrVar$dataTable.CV)
  expect_type(PLTdensOfVar,"list")


  PLTdisHIST=PlotDiscrHist(sampleA,"LDL")
  expect_s3_class(PLTdisHIST$layers[[1]], "ggproto")

  PLTcorrWL=PlotCorrWithRegrLine(sampleA,"CHOL", "HDL")
  expect_s3_class(PLTcorrWL$layers[[1]], "ggproto")

  LDL_empirVrnc = var(sampleA$LDL)
  LDL_errPropVrnc = LDLErrPrp(sampleA$CHOL,sampleA$HDL,sampleA$TG)
  LDLbootStrp=as.data.frame(LDLbootVrnc(sampleA$CHOL,sampleA$HDL,sampleA$TG,noOfReps = 2))
  LDLdensBoot=LDL_DensityPlotOfbootst(LDLbootStrp,"Title",LDL_empirVrnc,LDL_errPropVrnc)
  expect_error(LDLdensBoot, NA)

  sampleA$AIP = AIPcalc(sampleA$TG,sampleA$HDL, SI=FALSE)
  AIP_empirVrnc=var(sampleA$AIP)
  AIP_errPropVrnc=AIPErrPrp(sampleA$TG,sampleA$HDL, SI=FALSE)
  AIP_errPropVrnc2Ord=AIPErrPrp2Ord(sampleA$TG,sampleA$HDL, SI=FALSE)
  DfAIPboost=as.data.frame(AIPbootVrnc(sampleA$TG,sampleA$HDL,noOfReps = 2, SI=FALSE))
  AIPdensBoot=AIP_DensityPlotOfbootst(DfAIPboost,"Title",AIP_empirVrnc, AIP_errPropVrnc, AIP_errPropVrnc2Ord)
  expect_error(AIPdensBoot, NA)

})
