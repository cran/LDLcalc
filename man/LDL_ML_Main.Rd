% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/9ML_models_functions.R
\name{LDL_ML_Main}
\alias{LDL_ML_Main}
\title{Create, train, assess and return an ML prediction model.}
\usage{
LDL_ML_Main(DataCSV, partition, MLmethod, ReportMultiPlot = TRUE)
}
\arguments{
\item{DataCSV}{The .csv or a data table file, path containing the data with which the model will be trained and assessed. Must contain
at least 4 columns, named "CHOL", "HDL", "TG" and "LDLd", through which the train data and the validation data will be
extracted.}

\item{partition}{A value in the range (0,1) that stipulates what percentage of the input data will be
used for training the model, while the remainder will be used to assess it.}

\item{MLmethod}{A string that stipulates the Machine Learning method
("lm","rlm","glmnet","earth","svmRadial","knn","gbm","cubist" or "rf")
that is to be used to train the prediction model with.}

\item{ReportMultiPlot}{A boolean that allows the user to select whether the LDL_ML_Main function will
plot a diagram with 5 plots, relating different stats on the newly created model. Preset to TRUE.}
}
\value{
It initializes and returns the ML prediction model. In case of bad input,
it will return either -2 (illegitimate partition input) or -3 (illegitimate ML method input).
}
\description{
This function reads data from a DATACSV.csv, or a data table file. It partitions them according to the partition parameter
and labels them, trains the model (according to the ML method chosen and the first set of the partitioned data),
assesses the model using the second set of the partition data and returns it.
}
\examples{
model = LDL_ML_Main(SampleData,0.7,"lm",ReportMultiPlot=FALSE)
}
