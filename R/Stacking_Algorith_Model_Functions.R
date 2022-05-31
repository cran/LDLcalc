# This .R file requires "caret" package to run


# Using a Stacking Algorithm, creates a combination all models to predict LDL values
# This function creates and returns a prediction model made by stacking the 9 others in the LDL_ML_Main
#' @importFrom utils capture.output
#' @importFrom resample resample
LDL_ML_train_StackingAlgorithm<-function(trainData){
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
  StackingAlgorithm_Info<-summary(results)
  stackControl<-caret::trainControl(method="repeatedcv", number=10, index = caret::createFolds(trainData$LDLd, 10), repeats=3,
                             savePredictions=TRUE,search = "random")
  StackingAlgorithm_glm<-caretEnsemble::caretStack(models, method="glm", trControl=stackControl)
  return(list("stackModel" = StackingAlgorithm_glm, "results" = results))
}


#' Create, train, assess and return a Stacking Algorithm Machine Learning prediction model
#'
#' This function reads data from a DATACSV.csv or data table file. It partitions them according to the partition parameter
#' and labels them, trains all of the models and 'stacks' them into one, assesses them using the second set
#' of the partition data, optionally plots some info relating the accuracy of the models and returns them for further use.
#'
#' @param DataCSV The .csv or data table file, path containing the data with which the model will be trained and assessed. Must contain
#'  at least 4 columns, named "CHOL", "HDL", "TG" and "LDLd", through which the train data and the validation data will be
#'  extracted.
#' @param partition A value in the range (0,1) that stipulates what percentage of the input data will be
#'  used for training the model, while the remainder will be used to assess it.
#' @param ReportMultiPlot A boolean that allows the user to select whether the LDL_ML_Main function will
#'  plot a diagram with 5 plots, relating different stats on the newly created model. Preset to TRUE.
#' @param ComparisonPlot A boolean that allows the user to select whether the LDL_ML_Main_All_Models function will
#'  plot a comparison plot, relating different stats on the newly created models. Preset to TRUE.
#' @return It initializes and returns the stacked algorithm prediction model.
#' In case of bad input, it will return -2 (illegitimate partition input)
#' @examples
#' \donttest{
#' stackModel = LDL_ML_Main_StackingAlgorithm(SampleData,0.8,ReportMultiPlot=TRUE,ComparisonPlot=TRUE)
#' }
#' @export
LDL_ML_Main_StackingAlgorithm<-function(DataCSV, partition, ReportMultiPlot =TRUE, ComparisonPlot = TRUE){
  Data = CsvRead(DataCSV)
  if (partition<=0 || partition>=1){
    message("The partition is out of bounds. The limit are (0,1),and you entered : ")
    message (partition)
    return(-2)
  }

  result = dataPartitions(Data,partition)
  trainData = result$trainData
  validationData = result$validationData

  returned = LDL_ML_train_StackingAlgorithm(trainData)
  stackAlgModel = returned$stackModel
  results = returned$results

  stackAlgPredictions = LDL_ML_predict(stackAlgModel,validationData)

  if (ReportMultiPlot == TRUE){
    ReportMultiPlot(stackAlgPredictions,validationData$LDLd, "Stacking Algorithm Model")
  }
  if (ComparisonPlot == TRUE){
    a = Comparison_Models_Plot(results)
    print(a)

  }
  return(stackAlgModel)
}
