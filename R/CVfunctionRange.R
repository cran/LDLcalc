#' @title Change variance of distribution while keeping the mean constant
#' @description Takes as input a vector of values that have a distribution with a variance
#' expressed as coefficient of variation (CV), It then modifies the
#' distribution, so that the CV changes from a lower bound (lower_CV_bound)
#' to an upper bound (upper_CV_bound), while the mean stays constant.
#' Optionally, it can also plot the
#' resulting variances.
#' @param sampleVector The vector containing the initial values with initial
#' distribution and CV
#' @param lower_CV_Bound The lower CV value we wish the distribution to achieve.
#' @param upper_CV_Bound The upper CV value we wish the distribution to achieve.
#' @param maxRandIter The function uses a stochastic algorithm to change the
#' CV values from the lower to the upper value desired. If these bounds have not
#' been achieved after maxRandIter iterations (default=10000), the function will
#' exit, in order to avoid a possible infinite loop.
#' @param plot If set to TRUE, will plot the resulting variances (default=FALSE)
#' @return A data frame where each column contains a distribution of values with
#' increasing CV, from the lower bound to upper bound.
#' @examples
#' \dontrun{
#' DataFrame=CV_Range(sampleA$LDL,0,100,maxRandIter = 100, plot=TRUE)
#' }
#' @export
CV_Range <- function (sampleVector, lower_CV_Bound, upper_CV_Bound, maxRandIter=10000, plot=F){
  return(fChngDistrCVMeanConst(vecDistr = sampleVector, lowerCV = lower_CV_Bound, upperCV = upper_CV_Bound, maxIter=maxRandIter, plot=plot))
}


#' @title Calculate the CV of a set of values.
#' @description This function calculates the coefficient of variation (CV) of the values supplied.
#' @param Vec The vector of values for which to calculate the CV.
#' @param roundDigits Number of digits to round the result to (default=2).
#' @return Returns the value of the CV.
#' @importFrom stats sd
#' @examples
#' CV = CV(sampleA$LDL)
#' @export
CV <- function(Vec, roundDigits=2) {
  return(round((sd(Vec) /mean(Vec))*100, roundDigits))
}
