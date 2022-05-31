#' @title Plots LDL variance versus increasing cholesterol variance
#' @description Plots the variance of LDL (both Error Propagation and Bootstrap variance)
#' versus increasing cholesterol variance.
#' @param dfCholVars A data frame with increasing variances of cholesterol.
#' @param LDLCHOLVrncErrProp A vector with the error propagation variances of LDL when the
#' cholesterol variance increases.
#' @param LDLCHOLVrncBoot A vector with the bootstrap variances of LDL when the
#' cholesterol variance increases.
#' @return The function creates a plot with the LDL error propagation and
#' bootstrap variances versus increasing cholesterol variances.
#' @examples
#' \donttest{
#' # For HDL - CHOL:
#' # Make the DF of ascending variances for CHOL of sample
#' CHOLVariances = CV_Range(sampleA$CHOL,1,10,plot=FALSE)
#' LDLCHOLDependance = LDL_CHOLVrnc(CHOLVariances, sampleA$HDL, sampleA$TG, bootStrpReps=2000)
#' plotCholVrncToLDL(CHOLVariances,LDLCHOLDependance$ErrPropVrnc,LDLCHOLDependance$BootVrnc)
#' }
#' @export
plotCholVrncToLDL <- function(dfCholVars,LDLCHOLVrncErrProp,LDLCHOLVrncBoot){
  par(xpd=NA); # Sets plot clipping to the device.
  plot(apply(dfCholVars, 2, var), LDLCHOLVrncErrProp, type="l", lty=2,
       ylim=c(0, max(max(LDLCHOLVrncBoot), max(LDLCHOLVrncErrProp)) ),
       lwd=1.5, xlab="CHOL Variance", ylab="LDL Variance")
  title(paste("LDL variance relative to CHOL variance"), line=0.5)
  lines(apply(dfCholVars, 2, var), LDLCHOLVrncBoot, type="l", lty=1, lwd=1.5)
  legend("bottomright",  lty=c(2,1), legend=c("LDL Error Propagation Variance", "LDL Bootstrap Variance"))
}


#' @title Plots LDL variance versus increasing HDL variance
#' @description Plots the variance of LDL (both Error Propagation and Bootstrap variance)
#' versus increasing HDL variance.
#' @param dfHDLVars A data frame with increasing variances of HDL.
#' @param LDLHDLVrncErrProp A vector with the error propagation variances of LDL when the
#' HDL variance increases.
#' @param LDLHDLVrncBoot A vector with the bootstrap variances of LDL when the
#' HDL variance increases.
#' @return The function creates a plot with the LDL error propagation and
#' bootstrap variances versus increasing HDL variances.
#' @examples
#' \donttest{
#' ## For LDL - HDL:
#' # Make the DF of ascending variances for HDL of sample
#' HDLVariances = CV_Range(sampleA$HDL,15,25,plot=FALSE)
#' # Get the Error Propagation and the Bootstrap variance of LDL relative to HDL
#' LDLHDLSampleDependance = LDL_HDLVrnc(HDLVariances,sampleA$CHOL, sampleA$TG, bootStrpReps=2000)
#' plotHDLVrncToLDL(HDLVariances,LDLHDLSampleDependance$ErrPropVrnc,LDLHDLSampleDependance$BootVrnc)
#' }
#' @export
plotHDLVrncToLDL <- function(dfHDLVars,LDLHDLVrncErrProp, LDLHDLVrncBoot){
  par(xpd=NA); # Sets plot clipping to the device.
  plot(apply(dfHDLVars, 2, var), LDLHDLVrncErrProp, type="l", lty=2,
       ylim=c(0, max(max(LDLHDLVrncBoot), max(LDLHDLVrncErrProp)) ),
       lwd=1.5, xlab="HDL Variance", ylab="LDL Variance")
  title(paste("LDL variance relative to HDL variance"), line=0.5)
  lines(apply(dfHDLVars, 2, var), LDLHDLVrncBoot, type="l", lty=1, lwd=1.5)
  legend("bottomright",  lty=c(2,1), legend=c("LDL Error Propagation Variance", "LDL Bootstrap Variance"))
}

#' @title Plots LDL variance versus increasing triglyceride variance
#' @description Plots the variance of LDL (both Error Propagation and Bootstrap variance)
#' versus increasing triglyceride variance.
#' @param dfTGVars A data frame with increasing variances of triglycerides.
#' @param LDLTGVrncErrProp A vector with the error propagation variances of LDL when the
#' triglyceride variance increases.
#' @param LDLTGVrncBoot A vector with the bootstrap variances of LDL when the
#' triglycerides variance increases.
#' @return The function creates a plot with the LDL error propagation and
#' bootstrap variances versus increasing  triglyceride variances.
#' @examples
#' \donttest{
#' ## For LDL - TG:
#' # Make the DF of ascending variances for TG of sample
#' TGVariances = CV_Range(sampleA$TG,15,16,plot=FALSE)
#' # Get the Error Propagation and the Bootstrap variance of LDL relative to TG
#' LDLTGSampleDependance = LDL_TGVrnc(TGVariances,sampleA$CHOL, sampleA$HDL, bootStrpReps =2000)
#' plotTGVrncToLDL(TGVariances,LDLTGSampleDependance$ErrPropVrnc,LDLTGSampleDependance$BootVrnc)
#' }
#' @importFrom stats var
#' @export
plotTGVrncToLDL <- function(dfTGVars,LDLTGVrncErrProp, LDLTGVrncBoot){
  par(xpd=NA); # Sets plot clipping to the device.
  plot(apply(dfTGVars, 2, var), LDLTGVrncErrProp, type="l", lty=2,
       ylim=c(0, max(max(LDLTGVrncBoot), max(LDLTGVrncErrProp)) ),
       lwd=1.5, xlab="TG Variance", ylab="LDL Variance")
  title(paste("LDL variance relative to TG variance"), line=0.5)
  lines(apply(dfTGVars, 2, var), LDLTGVrncBoot, type="l", lty=1, lwd=1.5)
  legend("bottomright",  lty=c(2,1), legend=c("LDL Error Propagation Variance", "LDL Bootstrap Variance"))
}

#' @title Plots AIP variance versus increasing triglyceride variance
#' @description Plots the variance of AIP (both Error Propagation first and
#' second order, as well as Bootstrap variance) versus increasing triglyceride variance.
#' @param TGVars A data frame with increasing variances of triglycerides.
#' @param AIPVrncErrProp A vector with the first order error propagation
#' variances of AIP when the triglyceride variance increases.
#' @param AIPVrncErrProp2Ord A vector with the second order error propagation
#' variances of AIP when the triglyceride variance increases.
#' @param AIPbootVrnc A vector with the bootstrap
#' variances of AIP when the triglyceride variance increases.
#' @return The function creates a plot with the AIP error propagation and
#' bootstrap variances versus increasing  triglyceride variance.
#' @examples
#' \donttest{
#' # AIP - TG Variance
#' TGVariances = CV_Range(sampleA$TG,15,16,plot=FALSE)
#' AIPVrncChngTGVrnc = AIP_TGVrnc(TGVariances,sampleA$HDL,bootStrpReps = 2000)
#' TGErrPropVrnc = AIPVrncChngTGVrnc$ErrPropVrnc
#' TGErrPropVrnc2Ord = AIPVrncChngTGVrnc$ErrPropVrnc2Ord
#' TGBootVrnc = AIPVrncChngTGVrnc$BootVrnc
#' plotAIP_TGVrnc(TGVariances,TGErrPropVrnc,TGErrPropVrnc2Ord,TGBootVrnc)
#' }
#' @export
plotAIP_TGVrnc <- function(TGVars,AIPVrncErrProp, AIPVrncErrProp2Ord, AIPbootVrnc){
  par(xpd=NA); # Sets plot clipping to the device.
  plot(apply(TGVars, 2, var), AIPVrncErrProp, type="l", lty=1,
       ylim=c(0, max(max(AIPbootVrnc), max(AIPVrncErrProp)) ),
       lwd=1.5, xlab="HDL Variance", ylab="AIP Variance")
  title(paste("AIP relative to TG variance"), line=0.5)
  lines(apply(TGVars, 2, var), AIPVrncErrProp2Ord, type="l", lty=2)
  lines(apply(TGVars, 2, var), AIPbootVrnc, type="l", lty=3, lwd=1.5)
  legend("bottomright", lty=c(1,2,3), legend=c("AIP Error Propagation Variance",
                                               "AIP 2nd Order Taylor AIP Error Propagation Variance",
                                               "AIP Bootstrap Variance"))
}


#' @title Plots AIP variance versus increasing HDL variance
#' @description Plots the variance of AIP (both Error Propagation first and
#' second order, as well as Bootstrap variance) versus increasing HDL variance.
#' @param dfHDLVars A data frame with increasing variances of HDL.
#' @param AIPVrncErrProp A vector with the first order error propagation
#' variances of AIP when the HDL variance increases.
#' @param AIPVrncErrProp2Ord A vector with the second order error propagation
#' variances of AIP when the HDL variance increases.
#' @param AIPbootVrnc A vector with the bootstrap
#' variances of AIP when the triglyceride variance increases.
#' @return The function creates a plot with the AIP error propagation and
#' bootstrap variances versus increasing  triglyceride variance.
#' @importFrom graphics title lines legend
#' @examples
#' \donttest{
## AIP - HDL Variance
#' HDLVariances = CV_Range(sampleA$HDL,15,25,plot=FALSE)
#' AIPVrncChngHDLVrnc = AIP_HDLVrnc(HDLVariances,sampleA$TG, bootStrpReps=2000)
#' HDLErrPropVrnc = AIPVrncChngHDLVrnc$ErrPropVrnc
#' HDLErrPropVrnc2Ord = AIPVrncChngHDLVrnc$ErrPropVrnc2Ord
#' HDLBootVrnc = AIPVrncChngHDLVrnc$BootVrnc
#' plotAIP_HDLVrnc(HDLVariances,HDLErrPropVrnc,HDLErrPropVrnc2Ord,HDLBootVrnc)
#' }
#' @export
plotAIP_HDLVrnc <- function(dfHDLVars,AIPVrncErrProp,AIPVrncErrProp2Ord,AIPbootVrnc){
  par(xpd=NA); # Sets plot clipping to the device.
  plot(apply(dfHDLVars, 2, var), AIPVrncErrProp, type="l", lty=1,
       ylim=c(0, max(max(AIPbootVrnc), max(AIPVrncErrProp)) ),
       lwd=1.5, xlab="HDL Variance", ylab="AIP Variance")
  title(paste("AIP relative to HDL variance"), line=0.5)
  lines(apply(dfHDLVars, 2, var), AIPVrncErrProp2Ord, type="l", lty=2)
  lines(apply(dfHDLVars, 2, var), AIPbootVrnc, type="l", lty=3, lwd=1.5)
  legend("bottomright", lty=c(1,2,3), legend=c("AIP Error Propagation Variance",
                                        "AIP 2nd Order Taylor AIP Error Propagation Variance",
                                        "AIP Bootstrap Variance"))
}
