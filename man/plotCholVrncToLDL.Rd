% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotsK.R
\name{plotCholVrncToLDL}
\alias{plotCholVrncToLDL}
\title{Plots LDL variance versus increasing cholesterol variance}
\usage{
plotCholVrncToLDL(dfCholVars, LDLCHOLVrncErrProp, LDLCHOLVrncBoot)
}
\arguments{
\item{dfCholVars}{A data frame with increasing variances of cholesterol.}

\item{LDLCHOLVrncErrProp}{A vector with the error propagation variances of LDL when the
cholesterol variance increases.}

\item{LDLCHOLVrncBoot}{A vector with the bootstrap variances of LDL when the
cholesterol variance increases.}
}
\value{
The function creates a plot with the LDL error propagation and
bootstrap variances versus increasing cholesterol variances.
}
\description{
Plots the variance of LDL (both Error Propagation and Bootstrap variance)
versus increasing cholesterol variance.
}
\examples{
\donttest{
# For HDL - CHOL:
# Make the DF of ascending variances for CHOL of sample
CHOLVariances = CV_Range(sampleA$CHOL,1,10,plot=FALSE)
LDLCHOLDependance = LDL_CHOLVrnc(CHOLVariances, sampleA$HDL, sampleA$TG, bootStrpReps=2000)
plotCholVrncToLDL(CHOLVariances,LDLCHOLDependance$ErrPropVrnc,LDLCHOLDependance$BootVrnc)
}
}
