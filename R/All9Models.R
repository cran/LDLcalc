# This .R file requires "caret" and "resample" packages to run


#' Create, train, assess and return all ML prediction models
#'
#' This function reads data from a DATACSV.csv or a data.table file.It partitions them according to the partition parameter
#' and labels them, trains all of the models, assesses them using the second set of the partition data, optionally
#' plots some info relating the accuracy of the models and returns them for further use.
#'
#' @param DataCSV The .csv or a data table file, path containing the data with which the model will be trained and assessed. Must contain
#'  at least 4 columns, named "CHOL", "HDL", "TG" and "LDLd", through which the train data and the validation data will be
#'  extracted.
#' @param partition A value in the range (0,1) that stipulates what percentage of the input data will be
#'  used for training the models, while the remainder will be used to assess them.
#' @param ReportMultiPlot A boolean that allows the user to select whether the LDL_ML_Main function will
#'  plot a diagram with 5 plots, relating different stats on the newly created model. Preset to TRUE.
#' @param ComparisonPlot A boolean that allows the user to select whether the LDL_ML_Main_All_Models function will
#'  plot a comparison plot, relating different stats on the newly created models. Preset to TRUE.
#' @return It initializes and returns all ML prediction models. In case of bad input,
#'  it will return  -2 (illegitimate partition input).
#' @examples
#' \donttest{
#' allModels = LDL_ML_Main_All_Models(SampleData,0.8,ReportMultiPlot = FALSE,ComparisonPlot=FALSE)
#' }
#' @importFrom utils capture.output
#' @export
LDL_ML_Main_All_Models<-function(DataCSV, partition, ReportMultiPlot = TRUE, ComparisonPlot = TRUE){
  Data = CsvRead(DataCSV)
  if (partition<=0 || partition>=1){
    message("The partition is out of bounds. The limit are (0,1),and you entered : ")
    message(partition)
    return(-2)
  }
  result = dataPartitions(Data,partition)
  trainData = result$trainData
  validationData = result$validationData
  Folds=0
  dimtrain=dim(trainData)
  if (dimtrain[1]>=200){
    Folds=8
  }else {Folds=2}

  control <- caret::trainControl(method="repeatedcv", number=10, index = caret::createFolds(trainData$LDLd, Folds), repeats=3, savePredictions='final', classProbs=FALSE)
  StackingAlgorithm_List<-c("lm","rlm","glmnet","earth","svmRadial","knn","gbm","cubist","rf")

  invisible(capture.output(models<-caretEnsemble::caretList(LDLd ~ CHOL + HDL + TG, data=trainData, trControl=control,
                                   methodList=StackingAlgorithm_List)))
  results<-caret::resamples(models)
  Info<-summary(results)
  allModelsPredictions = LDL_ML_predict(models,validationData)
  modelNameList = c("Linear Model","Robust Linear Model","Generalized Linear Model Net","Earth Model","Support Vector Machine Radial Model","k-Nearest Neighbors Algorithm Model","Gradient Boosting Model","Cubist Model","Random Forest Model")
  if (ReportMultiPlot == TRUE){
    for (i in 1:9){
    ReportMultiPlot(allModelsPredictions[,i],validationData$LDLd,modelNameList[i])
    }
  }
  if (ComparisonPlot == TRUE){
    a = Comparison_Models_Plot(results)
    print(a)
  }

  return (models)
}
