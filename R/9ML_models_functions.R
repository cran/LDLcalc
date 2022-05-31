# This .R file requires "caret" libraries to run



# This function trains your ML-based prediction model with the train data you provide it, based on the ML
# method you choose.
LDL_ML_train<-function(trainData,MLmethod){
  if (MLmethod=="lm"){
    model_lm<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData,metric="RMSE",method="lm",
                           preProc=c("center", "scale"),trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_lm)}

  else if (MLmethod=="rlm"){
    model_rlm<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="rlm",
                            metric="RMSE", preProc=c("center", "scale"), trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 3))
    return(model_rlm)}

  else if (MLmethod=="glmnet"){
    model_glmnet<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="glmnet",
                               metric="RMSE", preProc=c("center", "scale"), trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_glmnet)}

  else if (MLmethod=="earth"){
    model_earth<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="earth",
                              metric="RMSE", preProc=c("center", "scale"),tuneGrid=expand.grid(.degree=1, .nprune=10), trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_earth)}

  else if (MLmethod=="svmRadial"){
    model_svmRadial<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="svmRadial",
                                  metric="RMSE", preProc=c("BoxCox"), trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_svmRadial)}

  else if (MLmethod=="knn"){
    model_knn<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="knn",
                            metric="RMSE",tuneGrid = expand.grid(k=2), preProc=c("BoxCox"), trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_knn)}

  else if (MLmethod=="gbm"){
    model_gbm<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="gbm",preProc=c("center", "scale"),
                            metric="RMSE",verbose=FALSE,trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_gbm)}

  else if (MLmethod=="cubist"){
    model_cubist<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="cubist",preProc=c("center", "scale"),
                               metric="RMSE", trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_cubist)}

  else if (MLmethod=="rf"){
    model_rf<-caret::train(LDLd ~ CHOL + HDL + TG, data=trainData, method="rf",preProc=c("center", "scale"),
                           metric="RMSE", trControl=caret::trainControl(method="repeatedcv",number=5,repeats = 5))
    return(model_rf)}
}

#' @title Predict LDL value(s)
#' @description This function predicts and returns  predictions, based on the model previously trained.
#' @param model The model with which the predictions will be made.
#' @param data The data with which the predictions will be made, can either be a single set of (CHOL,HDL,TG) values
#' or a data table of sets of said values.
#' @return The predicted LDL value(s).
#' @importFrom stats predict
#' @examples
#' modelPrediction = LDL_ML_predict(model$model,data.table::data.table(CHOL=170.5,HDL=35.12,TG=175))
#'
#' @export
LDL_ML_predict<-function(model,data){
  predictions <- predict(model,data)
  return(predictions)
}

#' @title  Create, train, assess and return an ML prediction model.
#' @description This function reads data from a DATACSV.csv, or a data table file. It partitions them according to the partition parameter
#' and labels them, trains the model (according to the ML method chosen and the first set of the partitioned data),
#' assesses the model using the second set of the partition data and returns it.
#' @param DataCSV The .csv or a data table file, path containing the data with which the model will be trained and assessed. Must contain
#' at least 4 columns, named "CHOL", "HDL", "TG" and "LDLd", through which the train data and the validation data will be
#' extracted.
#' @param partition A value in the range (0,1) that stipulates what percentage of the input data will be
#'  used for training the model, while the remainder will be used to assess it.
#' @param MLmethod A string that stipulates the Machine Learning method
#' ("lm","rlm","glmnet","earth","svmRadial","knn","gbm","cubist" or "rf")
#' that is to be used to train the prediction model with.
#' @param ReportMultiPlot A boolean that allows the user to select whether the LDL_ML_Main function will
#' plot a diagram with 5 plots, relating different stats on the newly created model. Preset to TRUE.
#' @return It initializes and returns the ML prediction model. In case of bad input,
#' it will return either -2 (illegitimate partition input) or -3 (illegitimate ML method input).
#' @examples
#' model = LDL_ML_Main(SampleData,0.7,"lm",ReportMultiPlot=FALSE)
#' @export
LDL_ML_Main<-function(DataCSV, partition, MLmethod, ReportMultiPlot = TRUE){
  Data = CsvRead(DataCSV)
  if (partition<=0 || partition>=1){
    message("The partition is out of bounds. The limit are (0,1),and you entered : ")
    message(partition)
    return(-2)
  }
  MLmethodList<-c("lm","rlm","glmnet","earth","svmRadial","knn","gbm","cubist","rf")
  if (!(MLmethod%in%MLmethodList)){
    message("ML method entered does not exist. You entered: ")
    message(MLmethod)
    return (-3)
  }

  result = dataPartitions(Data,partition)
  trainData = result$trainData
  validationData = result$validationData

  model = LDL_ML_train(trainData,MLmethod)
  predictions = LDL_ML_predict(model,validationData)

  if (ReportMultiPlot == TRUE){
    modelNamesList = c("Linear Model","Robust Linear Model","Generalized Linear Model Net","Earth Model","Support Vector Machine Radial Model","k-Nearest Neighbors Algorithm Model","Gradient Boosting Model","Cubist Model","Random Forest Model")
    for (i in 1:9)
      {
      if (MLmethod == MLmethodList[i]){
        index = i
        }
    }
    # If the selected ML method is "earth", the predictions vector will instead be returned as a 2D one column matrix
    # If we want to plot the predictions with ReportMultiPlot, we must first convert said matrix into a vector
    if (MLmethod == "earth"){
      predictions = as.vector(predictions)
    }
    ReportMultiPlot(predictions,validationData$LDLd,modelNamesList[index])
  }
  return(list("model"=model,"trainData"=trainData,"testData"=validationData))
}
