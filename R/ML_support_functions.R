# This .R file requires "data.table" and "caret" packages to run

# Read the .csv file
CsvRead<-function(CsvData){
  if (length(class(CsvData)) == 1 && class(CsvData) == "character"){
    Data<-data.table::fread(CsvData)
    return(Data)
  }
  else if (sum(ifelse(class(CsvData) == c("data.table","data.frame"),1,0)) == 2){
    return(CsvData)
  }
  else
  {
    stop("The data you entered is in neither .csv nor data.table format.")
  }
}

# Create partitions from your .csv file extracted data table
dataPartitions<-function(Data,partition){
  i<-caret::createDataPartition(Data$LDLd,p=partition,list=FALSE)
  trainData<-Data[i,]
  validationData <- Data[-i,]
  return(list("trainData" = trainData,"validationData" = validationData))
}

