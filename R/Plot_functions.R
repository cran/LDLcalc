# This .R file requires "lares", "corrplot" and "RColorBrewer" packages to run


# Creates a multi-plot diagram reporting very useful information about the model
ReportMultiPlot<-function(inputModelPredict,validationDataColumn, name = ""){
  lares::mplot_full(tag = validationDataColumn,
             score = inputModelPredict,
             splits =10 , subtitle = name)
}


# Create a plot to compare all models
Comparison_Models_Plot<-function(results){
  scales <- list(x=list(relation="free"), y=list(relation="free"))
  lattice::dotplot(results)
  lattice::bwplot(results, scales=scales)
  lattice::splom(results)
}


#' Plots a correlation plot to see the correlation between different columns of your data, for example LDL,HDL relative to age
#'
#' @param listDaten This is a list of the data amongst which one can see the correlation. One needs to provide at least two
#' columns of equal length in order to see the correlation between them
#' @importFrom stats cor
#' @return No return value, it prints the requested diagram.
#' @examples
#' Corr_Plot(SampleData[1:5])
#' @export
Corr_Plot<-function(listDaten){
  M <-cor(listDaten)
  corrplot::corrplot(M, type="upper", order="hclust",
           col=RColorBrewer::brewer.pal(n=8, name="RdYlBu"))
}
