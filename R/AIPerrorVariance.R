#' @title Calculate the variance of AIP when HDL changes.
#' @description This function calculates the variance of the Atherogenic Index of Plasma
#' (AIP) using two methods:
#' Error Propagation (first and second order) and Bootstrap when the HDL
#' distribution changes. It uses a
#' data frame where each column contains a different set of HDL values of
#' increasing variance as the column index increases.
#' @param TG A vector or data frame column containing the triglyceride (TG)
#' values to be used for the calculation of the variance of the
#' Atherogenic Index of Plasma (AIP).
#' @param dfHDL A data frame where each column contains a different set of HDL
#' values of increasing variance as the column index increases.
#' @param bootStrpReps (Default=2000) Number of bootstrap iterations
#' for bootstrap variance calculation.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI
#' units for TG and HDL (mmol/L). If mg/dl units are provided instead, SI must
#' be set to FALSE.
#' @return It returns a list with the first order Error Propagation variance
#' (ErrPropVrnc), the second order Error Propagation variance (ErrPropVrnc2Ord)
#' and the bootstrap variance (BootVrnc). Each list element is a vector of
#' length equal to the number of columns of the HDL data frame supplied as
#' argument (dfHDL) and each vector value corresponds to the respective
#' variance of the corresponding data frame columns.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' dfHDL = CV_Range(sampleB$HDL,0,100,maxRandIter = 1000, plot=FALSE)
#' AIP_HDLVrnc=AIP_HDLVrnc(dfHDL,sampleA$TG, bootStrpReps=2000)
#' }
#' @export
AIP_HDLVrnc <- function(dfHDL, TG, SI=TRUE, bootStrpReps=2000) {
  if (SI==FALSE){
    TG=TG*0.01129
  }
  AIPVrncChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfHDL))
  AIPVrnc2OrdChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfHDL))
  AIPVrncChangingBootVrnc = vector(mode="numeric", length=ncol(dfHDL))
  pb = txtProgressBar(min = 0, max = ncol(dfHDL), initial = 0, style=3)
  for(i in 1:ncol(dfHDL)) {
    AIPVrncErrProp = fAIPErrPrp(TG,  dfHDL[,i])
    AIPVrncErrProp2Ord = fAIPErrPrp2Ord(TG,  dfHDL[,i])
    AIPVrncChangingErrPropVrnc[i] = AIPVrncErrProp
    AIPVrnc2OrdChangingErrPropVrnc[i] = AIPVrncErrProp2Ord
    AIPVrncBoot = AIPbootVrnc(TG, dfHDL[,i], noOfReps = bootStrpReps)
    AIPVrncChangingBootVrnc[i] = AIPVrncBoot$Var
    setTxtProgressBar(pb, i)
  }
  return(list(ErrPropVrnc = AIPVrncChangingErrPropVrnc,
              ErrPropVrnc2Ord = AIPVrnc2OrdChangingErrPropVrnc,
              BootVrnc = AIPVrncChangingBootVrnc))
}


#' @title Calculate the variance of AIP when TG changes.
#' @description Calculate the variance of the Atherogenic Index of Plasma (AIP)
#' using two methods:
#' Error Propagation (first and second order) and Bootstrap when the
#' triglyceride distribution changes. It uses a
#' data frame where each column contains a different set of TG values of
#' increasing variance as the column index increases.
#' @param HDL A vector or data frame column containing the HDL
#' values to be used for the calculation of the variance of the
#' Atherogenic Index of Plasma (AIP).
#' @param dfTG A data frame where each column contains a different set of TG
#' values of increasing variance as the column index increases.
#' @param bootStrpReps (Default=2000) Number of bootstrap iterations
#' for bootstrap variance calculation.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI
#' units for TG and HDL (mmol/L). If mg/dl units are provided instead, SI must
#' be set to FALSE.
#' @return It returns a list with the first order Error Propagation variance
#' (ErrPropVrnc), the second order Error Propagation variance (ErrPropVrnc2Ord)
#' and the bootstrap variance (BootVrnc). Each list element is a vector of
#' length equal to the number of columns of the TG data frame supplied as
#' argument (dfTG) and each vector value corresponds to the respective
#' variance of the corresponding data frame columns.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' dfTG = CV_Range(sampleA$TG, 0, 100, maxRandIter = 100, plot=FALSE)
#' AIP_TGVrnc=AIP_TGVrnc(dfTG,sampleA$HDL, bootStrpReps=2000)
#' }
#' @export
AIP_TGVrnc <- function(dfTG, HDL, SI=TRUE, bootStrpReps=2000) {
  if (SI==FALSE){
    HDL=HDL*0.0259
  }
  AIPVrncChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfTG))
  AIPVrnc2OrdChangingErrPropVrnc = vector(mode="numeric", length=ncol(dfTG))
  AIPVrncChangingBootVrnc = vector(mode="numeric", length=ncol(dfTG))
  pb = txtProgressBar(min = 0, max = ncol(dfTG), initial = 0, style=3)
  for(i in 1:ncol(dfTG)) {
    AIPVrncErrProp = fAIPErrPrp(dfTG[,i], HDL)
    AIPVrncErrProp2Ord = fAIPErrPrp2Ord(dfTG[,i], HDL)
    AIPVrncChangingErrPropVrnc[i] = AIPVrncErrProp
    AIPVrnc2OrdChangingErrPropVrnc[i] = AIPVrncErrProp2Ord
    AIPVrncBoot = AIPbootVrnc(dfTG[,i], HDL, noOfReps = bootStrpReps)
    AIPVrncChangingBootVrnc[i] = AIPVrncBoot$Var
    setTxtProgressBar(pb, i)
  }
  return(list(ErrPropVrnc = AIPVrncChangingErrPropVrnc,
              ErrPropVrnc2Ord = AIPVrnc2OrdChangingErrPropVrnc,
              BootVrnc = AIPVrncChangingBootVrnc))
}


