# This .R file requires "data.table" and "caret" packages to run

# Read the .csv file
CsvRead<-function(CsvData){
  if (inherits(CsvData,"data.table")==TRUE){
    return(CsvData)
  }
  else if (inherits(CsvData,"data.table")==FALSE)
  {
    CsvData<-data.table::as.data.table(CsvData)
    return(CsvData)
  }
  else{
    stop("The data you entered is not a correct format.")
  }

}

# Create partitions from your .csv file extracted data table
dataPartitions<-function(Data,partition){
  i<-caret::createDataPartition(Data$LDLd,p=partition,list=FALSE)
  trainData<-Data[i,]
  validationData <- Data[-i,]
  return(list("trainData" = trainData,"validationData" = validationData))
}

