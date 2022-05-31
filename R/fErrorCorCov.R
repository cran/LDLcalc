#' @title Calculate correlation and covariance of residuals
#' @description This function calculates the residuals
#' (the squared difference of each value from the mean) of two groups. Then
#' it calculates the correlation and covariance between the residuals of the two
#' groups and plots them.
#' @param dataset1 A vector containing the values of the first group.
#' @param dataset2 A vector containing the values of the second group.
#' @param plot (default=TRUE). Plot the errors of the first group versus the second one.
#' @return Returns a list with the residuals of the two groups (distr1Error,
#' distr2Error), the correlation coefficient between the two groups (Correlation)
#' and the covariance of the two groups (Covariance).
#' @importFrom stats cov
#' @examples
#' \dontrun{
#' ErrorOFCorCov=ErrorCorCov(sampleA$HDL[1:20],sampleA$CHOL[1:20],plot = FALSE)
#'  }
#' @export
  ErrorCorCov <- function(dataset1,dataset2 , plot=F) {
  if(length(dataset1)!=length(dataset2)){message("The length of dataset1 and dataset2 the must be the same ")}
  else{
  dataset1Error = (dataset1 - mean(dataset1))^2
  dataset2Error = (dataset2 - mean(dataset2))^2
  correlation = round(cor(dataset1Error, dataset2),2)
  covariance = round(cov(dataset1Error, dataset2),2)}
  if(plot) {
    title = paste("Correlation between the errors of dataset1 and dataset2","\n Correlation Coefficient=",correlation, "\nCovariance=", covariance)
    plot(dataset1Error, dataset2Error, main=title, xlab=paste("dataset 1","Error"), ylab=paste("dataset 2", "Error"))
  }
  return(list(distr1Error=dataset1Error, distr2Error=dataset2Error, Correlation=correlation, Covariance=covariance))
}

