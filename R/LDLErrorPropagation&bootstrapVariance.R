# This .R file requires the "data.table" package to run.

#' @title Calculate variance of LDL using bootstrapping.
#' @description Function to calculate the variance of LDL using Bootstrapping.
#' @param CHOL A vector containing the cholesterol values to be used for LDL calculation.
#' @param HDL A vector containing the HDL values to be used for LDL calculation.
#' @param TG A vector containing the triglyceride values to be used for LDL calculation.
#' The three vectors (vecCHOL, vecHDL and vecG) must contain the same number of values
#' @param sampleSize Number of samples drawn uniformly and with replacement
#' from the three vectors. It cannot be larger than the sample size
#' (default=sample size).
#' @param noOfReps Number of bootstrap iterations.
#' @param pb Draw a text progress bar (defaut=FALSE)
#' @return It returns a data table with four columns. The first
#' column contains the mean of the LDL values for each iteration. The second
#' column contains the median of each iteration.The third
#' column contains the variance and the fourth column contains
#' the CV of each iteration. It also returns the median of the "Mean", "Var"
#' and "CV" columns of the data table.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' LDLboostrpVar = LDLbootVrnc(sampleA$CHOL, sampleA$HDL, sampleA$TG)
#' }
#' @export
  LDLbootVrnc <- function(CHOL, HDL, TG, sampleSize=length(CHOL), noOfReps=1000, pb=F) {
  if((length(CHOL) != length(HDL)) |
     (length(CHOL) != length(TG)) |
     (length(HDL) != length(TG)) ) {
    print("The three parameters must have the same number of measurements");
    return()
  }
  n = length(CHOL)
  if(sampleSize>n) {
    print("Size of bootstraped datasets cannot be larger than the original. Setting k=n")
    sampleSize=n;
    print(paste("=", n))
  }
  dfTmp <- data.frame(cbind(CHOL, HDL, TG))
  colnames(dfTmp) <- c("CHOL", "HDL", "TG")
  dtLDLBoot = data.table::data.table("Mean"=numeric(noOfReps),
                         "Median"=numeric(noOfReps),
                         "Var"=numeric(noOfReps),
                         "CV"=numeric(noOfReps))
  if(pb) pbar <- txtProgressBar(min = 0, max = noOfReps, style = 3)
  for(bootIdx in seq(1,noOfReps)) {
    dfSmpl <- dfTmp[sample(nrow(dfTmp), size=sampleSize, replace = T),]
    dfLDLSmpl = dfSmpl$CHOL - dfSmpl$HDL - (dfSmpl$TG/5)
    LDLSmplMean = mean(dfLDLSmpl)
    LDLSmplMedian = median(dfLDLSmpl)
    LDLSmplVar = var(dfLDLSmpl)
    LDLSmplCV = CV(dfLDLSmpl)
    sets::set(dtLDLBoot, i=bootIdx, j=1L, value=LDLSmplMean)
    sets::set(dtLDLBoot, i=bootIdx, j=2L, value=LDLSmplMedian)
    sets::set(dtLDLBoot, i=bootIdx, j=3L, value=LDLSmplVar)
    sets::set(dtLDLBoot, i=bootIdx, j=4L, value=LDLSmplCV)
    if(pb) setTxtProgressBar(pbar, bootIdx)
  }
  return(list("dataTable"=dtLDLBoot, "Mean"=median(dtLDLBoot$Mean),
              "Var"=median(dtLDLBoot$Var), "CV"=median(dtLDLBoot$CV)))
}

#' @title Calculate LDL variance using error propagation
#' @description Function to calculate the LDL Variance according to error propagation (delta) method.
#' @param CHOL A vector containing the cholesterol values to be used for LDL calculation.
#' @param HDL A vector containing the HDL values to be used for LDL calculation.
#' @param TG A vector containing the triglyceride values to be used for LDL calculation.
#' @param divFactor The factor by which to divide the triglyceride values
#' so as to approximate the VLDL. Default is 5, according to the Friedewald
#' equation.
#' @return The function returns the error propagation variance of LDL calculated
#' from the cholesterol, HDL and tryglyceride values passed as arguments.
#' @references  Casella G, Berger RL. Statistical Inference. 2nd ed. Duxbury Thomson Learning; 2002.
#' @importFrom stats cov
#' @examples
#' \dontrun{
#' LDLerrorPrp = LDLErrPrp(sampleA$CHOL, sampleA$HDL, sampleA$TG)
#' }
#' @export
LDLErrPrp <- function(CHOL, HDL, TG, divFactor=5) {
  LDLVar = var(CHOL) + var(HDL) + (var(TG)/(divFactor^2)) -2*cov(CHOL, HDL) - 2*(cov(CHOL,TG)/divFactor) + 2 * (cov(HDL,TG)/divFactor)
  return(LDLVar)
}

#' @title Calculate LDL variance when cholesterol variance changes
#' @description This function calculates the variance of LDL using two methods:
#' Error Propagation and Bootstrap when the cholesterol distribution changes. It uses a
#' data frame where each column contains a different set of cholesterol values of
#' increasing variance as the column index increases.
#' @param dfCHOL A data frame where each column contains a different set of
#' cholesterol values of increasing variance as the column index increases.
#' @param HDL A vector or data frame column containing the HDL
#' values to be used for the calculation of the variance of LDL.
#' @param TG A vector or data frame column containing the triglyceride
#' values to be used for the calculation of the variance of LDL.
#' @param bootStrpReps (Default=2000) Number of bootstrap iterations
#' for bootstrap variance calculation.
#' @return The function returns a list with the Error Propagation variance
#' (ErrPropVrnc) and the bootstrap variance (BootVrnc).
#' Each list element is a vector of
#' length equal to the number of columns of the cholesterol data frame supplied as
#' argument (dfCHOL) and each vector value corresponds to the respective
#' variance of the corresponding data frame columns.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' dfCHOL = CV_Range(sampleA$CHOL, 0, 100, maxRandIter = 3000, plot=FALSE)
#' LDLCHOLVar = LDL_CHOLVrnc(dfCHOL, sampleA$HDL, sampleA$TG, bootStrpReps=200)
#' }
#' @export
  LDL_CHOLVrnc <- function(dfCHOL, HDL, TG, bootStrpReps=2000) {
  LDLVrncChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfCHOL))
  LDLVrncChangingBootVrnc = vector(mode="numeric", length=ncol(dfCHOL))
  pb = txtProgressBar(min = 0, max = ncol(dfCHOL), initial = 0, style=3)
  for(i in 1:ncol(dfCHOL)) {
    LDLVrncErrProp = LDLErrPrp(dfCHOL[,i], HDL, TG)
    LDLVrncChangingErrPropVrnc[i] = LDLVrncErrProp
    LDLVrncBoot = LDLbootVrnc(dfCHOL[,i], HDL, TG, noOfReps = bootStrpReps)
    LDLVrncChangingBootVrnc[i] = LDLVrncBoot$Var
    setTxtProgressBar(pb, i)
  }
  return(list(ErrPropVrnc = LDLVrncChangingErrPropVrnc, BootVrnc = LDLVrncChangingBootVrnc))
}


#' @title Calculate LDL variance when HDL variance changes
#' @description This function calculates the variance of LDL using two methods:
#' Error Propagation and Bootstrap when the HDL distribution changes. It uses a
#' data frame where each column contains a different set of HDL values of
#' increasing variance as the column index increases.
#' @param dfHDL A data frame where each column contains a different set of
#' HDL values of increasing variance as the column index increases.
#' @param CHOL A vector or data frame column containing the cholesterol
#' values to be used for the calculation of the variance of LDL.
#' @param TG A vector or data frame column containing the triglyceride
#' values to be used for the calculation of the variance of LDL.
#' @param bootStrpReps (Default=2000) Number of bootstrap iterations
#' for bootstrap variance calculation.
#' @return The function returns a list with the Error Propagation variance
#' (ErrPropVrnc) and the bootstrap variance (BootVrnc).
#' Each list element is a vector of
#' length equal to the number of columns of the HDL data frame supplied as
#' argument (dfHDL) and each vector value corresponds to the respective
#' variance of the corresponding data frame columns.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' dfHDL = CV_Range(sampleA$HDL, 0, 100, maxRandIter = 4000, plot=FALSE)
#' LDLHDLVar = LDL_HDLVrnc(dfHDL,sampleA$CHOL,sampleA$TG,bootStrpReps=2000)
#' }
#' @export
LDL_HDLVrnc <- function(dfHDL, CHOL, TG, bootStrpReps=2000) {
  LDLVrncChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfHDL))
  LDLVrncChangingBootVrnc = vector(mode="numeric", length=ncol(dfHDL))
  pb = txtProgressBar(min = 0, max = ncol(dfHDL), initial = 0, style=3)
  for(i in 1:ncol(dfHDL)) {
    LDLVrncErrProp = LDLErrPrp(CHOL, dfHDL[,i], TG)
    LDLVrncChangingErrPropVrnc[i] = LDLVrncErrProp
    LDLVrncBoot = LDLbootVrnc(CHOL, dfHDL[,i], TG, noOfReps = bootStrpReps)
    LDLVrncChangingBootVrnc[i] = LDLVrncBoot$Var
    setTxtProgressBar(pb, i)
  }
  return(list(ErrPropVrnc = LDLVrncChangingErrPropVrnc, BootVrnc = LDLVrncChangingBootVrnc))
}


#' @title Calculate LDL variance when triglyceride variance changes
#' @description This function calculates the variance of LDL using two methods:
#' Error Propagation and Bootstrap when the triglyceride distribution changes.
#' It uses a data frame where each column contains a different set of triglyceride values of
#' increasing variance as the column index increases.
#' @param CHOL A vector or data frame column containing the cholesterol
#' values to be used for the calculation of the variance of LDL.
#' @param HDL A vector or data frame column containing the HDL
#' values to be used for the calculation of the variance of LDL.
#' @param dfTG A data frame where each column contains a different set of
#' triglyceride values of increasing variance as the column index increases.
#' @param bootStrpReps (Default=2000) Number of bootstrap iterations
#' for bootstrap variance calculation.
#' @return The function returns a list with the Error Propagation variance
#' (ErrPropVrnc) and the bootstrap variance (BootVrnc).
#' Each list element is a vector of
#' length equal to the number of columns of the triglyceride data frame supplied as
#' argument (dfTG) and each vector value corresponds to the respective
#' variance of the corresponding data frame columns.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' dfTG = CV_Range(sampleA$TG, 0, 100, maxRandIter = 3000, plot=FALSE)
#' LDLTGVar=LDL_TGVrnc( sampleA$CHOL,  sampleA$HDL, dfTG, bootStrpReps=2000)
#' }
#' @export
LDL_TGVrnc <- function(dfTG, CHOL, HDL, bootStrpReps=2000) {
  LDLVrncChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfTG))
  LDLVrncChangingBootVrnc = vector(mode="numeric", length=ncol(dfTG))
  pb = txtProgressBar(min = 0, max = ncol(dfTG), initial = 0, style=3)
  for(i in 1:ncol(dfTG)) {
    LDLVrncErrProp = LDLErrPrp(CHOL, HDL, dfTG[,i])
    LDLVrncChangingErrPropVrnc[i] = LDLVrncErrProp
    LDLVrncBoot = LDLbootVrnc(CHOL, HDL, dfTG[,i], noOfReps = bootStrpReps)
    LDLVrncChangingBootVrnc[i] = LDLVrncBoot$Var
    setTxtProgressBar(pb, i)
  }
  return(list(ErrPropVrnc = LDLVrncChangingErrPropVrnc, BootVrnc = LDLVrncChangingBootVrnc))
}


