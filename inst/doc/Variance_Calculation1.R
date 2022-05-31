## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
setwd(".")
require(LDLcalc)
require(ggplot2)
require(gridExtra)
require(data.table)
require(tidyr)

## ---- loadSample--------------------------------------------------------------

dfSmplA <- sampleA

## ---- calcAIP-----------------------------------------------------------------
dfSmplA$AIPrep1 <- AIPcalc(TG=dfSmplA$TGrep1,
                                    HDL=dfSmplA$HDLrep1, SI=F)
dfSmplA$AIPrep2 <- AIPcalc(TG=dfSmplA$TGrep2,
                                    HDL=dfSmplA$HDLrep2, SI=F)
dfSmplA$AIP <- AIPcalc(TG=dfSmplA$TG,
                                HDL=dfSmplA$HDL, SI=F)

## ---- discrHist, fig.height=4, fig.width=8, fig.align="center", fig.cap="Discrete histogram and the corresponding normal distribution."----
# Create a list of the parameters to be plotted.
lstParam <- list("CHOL", "HDL", "TG")
# Use lapply to apply the PlotDiscrHist function to all parameters.
lstPlt <- lapply(lstParam, LDLcalc::PlotDiscrHist, DF=dfSmplA)
do.call("grid.arrange", c(lstPlt, ncol = 3))

## ---- echo=F------------------------------------------------------------------
rm(lstParam, lstPlt)

## ---- fig.cap="Discrete histogram and corresponding normal distribution of LDL."----
PlotDiscrHist(dfSmplA, param="LDL", title="LDL")
JSDNormal(dfSmplA, "LDL")$JSD

## -----------------------------------------------------------------------------
distr1 <- seq(from=100, to=130, by=1)
distr2 <- seq(from=140, to=170, by=1)
JSD(distr1, distr2)$JSD

## ---- echo=F------------------------------------------------------------------
rm(distr1, distr2)

## ---- LDLchebyshev------------------------------------------------------------
LDLchebyshev <- chebyshev(vec=dfSmplA$LDL)
LDLchebyshev <- unlist(LDLchebyshev)
print(LDLchebyshev)

## ---- chebysevRange-----------------------------------------------------------
sum(dfSmplA$LDL>=157 & dfSmplA$LDL<=198) / nrow(dfSmplA) * 100

## ---- corcov------------------------------------------------------------------
lstParam1 <- list("CHOL", "CHOL", "HDL")# List of the first parameter to calculate correlation on.
lstParam2 <- list("HDL", "TG", "TG")# List of the second parameter to calculate correlation on.
#Function to calculate correlation and covariance of the pairs of parameters
fCor <- function(param1, param2) {
  correlation <- round(cor(dfSmplA[, param1], dfSmplA[, param2]), 2)
  covariance <- round(cov(dfSmplA[, param1], dfSmplA[, param2]), 2)
  return(list(correlation=correlation, covariance=covariance))
}

mtrxCorCov <- mapply(fCor, lstParam1, lstParam2)
colnames(mtrxCorCov) <- c("CHOL-HDL", "CHOL-TG", "HDL-TG")

print(mtrxCorCov)

## ---- coorplot, fig.height=4, fig.width=8, fig.align="center", fig.cap="Scatterplot of pairs of the parameters."----
grid.arrange(PlotCorrWithRegrLine(dfSmplA, "CHOL", "HDL"),
             PlotCorrWithRegrLine(dfSmplA, "CHOL", "TG"),
             PlotCorrWithRegrLine(dfSmplA, "HDL", "TG"),
             ncol=3)

## ---- empirVar----------------------------------------------------------------
vecEmpirVar <- c("LDL Empirical Variance" = round(var(dfSmplA$LDL), 2),
                 "AIP Empirical Variance" = round(var(dfSmplA$AIP), 4))
print(vecEmpirVar)

## ---- LDLErrPrpVar------------------------------------------------------------
LDLErrPrpVar <- LDLErrPrp(CHOL = dfSmplA$CHOL,
                                    HDL = dfSmplA$HDL,
                                    TG = dfSmplA$TG)
print(round(LDLErrPrpVar), 2)

## ---- AIPErrPrpVar------------------------------------------------------------
AIPErrPrpVar1Ord <-
  AIPErrPrp(TG = dfSmplA$TG,
                       HDL = dfSmplA$HDL,
                       SI = F)

AIPErrPrpVar2Ord <-
  AIPErrPrp2Ord(TG = dfSmplA$TG,
                       HDL = dfSmplA$HDL,
                       SI = F)

AIPErrProp <- c("AIP First Order Error Propagation" = AIPErrPrpVar1Ord,
                "AIP Second Order Error Propagation" = AIPErrPrpVar2Ord)
print(AIPErrProp)

## ---- echo=F------------------------------------------------------------------
vecVariances <- c(NA,
                  format(var(dfSmplA$LDL), nsmall=2),
                  format(LDLErrPrpVar, nsmall=2),
                  NA,
                  format(round(var(dfSmplA$AIP),5), nsmall=5),
                  format(AIPErrPrpVar1Ord, nsmall=5),
                  format(AIPErrPrpVar2Ord, nsmall=5))
dfVariances <- 
  cbind.data.frame(c("LDL",
                    "Empirical Variance",
                    "Error Propagation Variance",
                    "AIP",
                    "Empirical Variance",
                    "First-Order Error Propagation",
                    "Second-Order Error Propagation"),
    vecVariances)
colnames(dfVariances) <- c("Variance type", "Variance Value")

## ---- echo=F------------------------------------------------------------------
options(knitr.kable.NA = "")
knitr::kable(dfVariances, format = "latex", longtable = F, 
             caption="Empirical and Error Propagation Variance for LDL and AIP")

## ---- cache=F-----------------------------------------------------------------
LDLbootVar <- LDLbootVrnc(CHOL = dfSmplA$CHOL, HDL = dfSmplA$HDL,
                       TG = dfSmplA$TG, noOfReps = 10, pb=F)

AIPbootVar <- AIPbootVrnc(TG = dfSmplA$TG, HDL = dfSmplA$HDL,
                           SI= F, noOfReps = 10, pb=F)

## ---- echo=F------------------------------------------------------------------
vecVariances <- c(NA,
                  format(var(dfSmplA$LDL), nsmall=2),
                  format(LDLErrPrpVar, nsmall=2),
                  format(median(LDLbootVar$Var), nsmall=2),
                  NA,
                  format(round(var(dfSmplA$AIP),5), nsmall=5),
                  format(AIPErrPrpVar1Ord, nsmall=5),
                  format(AIPErrPrpVar2Ord, nsmall=5),
                  format(round(median(AIPbootVar$Var),5), nsmall=5))

dfVariances <- 
  cbind.data.frame(c("LDL",
                    "Empirical Variance",
                    "Error Propagation Variance",
                    "Bootstrap Variance",
                    "AIP",
                    "Empirical Variance",
                    "First-Order Error Propagation",
                    "Second-Order Error Propagation",
                    "Bootstrap Variance"),
    vecVariances)
colnames(dfVariances) <- c("Variance type", "Variance Value")

## ---- echo=F------------------------------------------------------------------
options(knitr.kable.NA = "")
knitr::kable(dfVariances, format = "latex", longtable = F,
             caption="Empirical, Error Propagation and Bootstrap Variance
             for LDL and AIP")

## ----eval=FALSE,fig.cap="Distribution density of the bootstrap LDL variance. The median is shown as a continuous vertical line and the 2.5 and 97.5 percentiles as dashed vertical lines. The empirical variance is depicted as a dotted vertical segment from the bottom to the middle of the graph. The error propagation variance is depicted as a dash-dot vertical segment from the top to the middle."----
#  LDL_DensityPlotOfbootst(LDLbootVar$dataTable,
#                                   empirVrnc = vecEmpirVar[[1]],
#                                   errPropVrnc = LDLErrPrpVar)

## ----eval=FALSE, fig.cap="Distribution density of the bootstrap AIP variance. The median is shown as a continuous vertical line and the 2.5 and 97.5 percentiles as dashed vertical lines. The empirical variance is depicted as a dotted vertical segment from the bottom to the middle of the graph. The error propagation variance is depicted as a dash-dot vertical segment from the top to the middle and the second order error propagation variance as a dashed segment in the middle."----
#  AIP_DensityPlotOfbootst(AIPbootVar$dataTable,
#                                   empirVrnc = vecEmpirVar[[2]],
#                                   errPropVrnc = AIPErrPrpVar1Ord,
#                                   errPropVrnc2Ord = AIPErrPrpVar2Ord)

## -----------------------------------------------------------------------------
CV(dfSmplA$CHOL)

## ---- CHOLVrncChng, fig.height=6,,fig.width=8, fig.align="center", position="!h", fig.cap="Means and variances for the CHOL distributions, in order to check that the mean stays constant and the variance increases.", cache=F----
dfSmplACHOLChngCV <- CV_Range(dfSmplA$CHOL,
                                       lower_CV_Bound = 0,
                                       upper_CV_Bound = 20,
                                       maxRandIter = 1000,
                                       plot=T)

## ---- LDL_CHOLVrnc500,results='hide',eval=FALSE, cache=F----------------------
#  LDL_Vrnc_Chng_CHOL <- LDL_CHOLVrnc(dfSmplACHOLChngCV,
#                                              dfSmplA$HDL,
#                                              dfSmplA$TG,
#                                              bootStrpReps = 500)

## ---- fig.height=3, fig.width=8,eval=FALSE, fig.cap="Variance of LDL (error propagation and boostrap for 500 iterations), when the CHOL variance changes."----
#  # First we will convert the list output of LDLcalc::LDL_CHOLVrnc to a
#  # data frame
#  LDL_Vrnc_Chng_CHOL <- do.call(cbind.data.frame, LDL_Vrnc_Chng_CHOL)
#  
#  # Create a data frame with the changing CHOL varince and the error propagation
#  # and bootstrap variances as columns, in order to use it for ggplot and
#  # remove the original data frame.
#  dfCHOLChngVrnc <- cbind.data.frame(CHOLVrnc= apply(dfSmplACHOLChngCV, 2, var),
#                   LDLErrPrp = LDL_Vrnc_Chng_CHOL$ErrPropVrnc,
#                   LDLBootstrp = LDL_Vrnc_Chng_CHOL$BootVrnc)
#  rm(LDL_Vrnc_Chng_CHOL)
#  
#  dfCHOLChngVrnc <- pivot_longer(dfCHOLChngVrnc, cols=2:3,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  
#  ggplot(data=dfCHOLChngVrnc,
#         aes(x=CHOLVrnc, y=Variance, linetype=Variance_Type)) +
#    geom_line() +
#    xlab("Cholesterol Variance") + ylab("LDL Variance") +
#    scale_linetype_discrete(name="Variance Type",
#                            labels = c("Bootstrap Variance",
#                                       "Error Propagation Variance"))

## ---- LDL_CHOLVrnc2000, results='hide',eval=FALSE,cache=F---------------------
#  LDL_Vrnc_Chng_CHOL <- LDL_CHOLVrnc(dfSmplACHOLChngCV,
#                                              dfSmplA$HDL,
#                                              dfSmplA$TG,
#                                              bootStrpReps = 2000)

## ---- echo=T, fig.height=3, fig.width=8,eval=FALSE, fig.cap="Variance of LDL (error propagation and boostrap for 2000 iterations), when the CHOL variance changes."----
#  LDL_Vrnc_Chng_CHOL <- do.call(cbind.data.frame, LDL_Vrnc_Chng_CHOL)
#  
#  dfCHOLChngVrnc <- cbind.data.frame(CHOLVrnc= apply(dfSmplACHOLChngCV, 2, var),
#                   LDLErrPrp = LDL_Vrnc_Chng_CHOL$ErrPropVrnc,
#                   LDLBootstrp = LDL_Vrnc_Chng_CHOL$BootVrnc)
#  rm(LDL_Vrnc_Chng_CHOL)
#  
#  dfCHOLChngVrnc <- pivot_longer(dfCHOLChngVrnc, cols=2:3,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  pltLDL_CHOLVrnc <-
#    ggplot(data=dfCHOLChngVrnc,
#           aes(x=CHOLVrnc, y=Variance, linetype=Variance_Type)) +
#      geom_line() +
#      xlab("Cholesterol Variance") + ylab("LDL Variance") +
#      scale_linetype_discrete(name="Variance Type",
#                              labels = c("Bootstrap Variance",
#                                         "Error Propagation Variance"))
#  pltLDL_CHOLVrnc

## -----------------------------------------------------------------------------
CV(dfSmplA$HDL)

## ---- HDLVrncChng, fig.height=7, fig.width=8, fig.align="center", fig.cap="Means and variances for the HDL distributions, in order to check that the mean stays constant and the variance increases.", cache=F----
dfSmplAHDLChngCV <- CV_Range(dfSmplA$HDL,
                                       lower_CV_Bound = 0,
                                       upper_CV_Bound = 20,
                                       maxRandIter = 10000,
                                       plot=T)

## ---- LDL_HDLVrnc2000,results='hide',eval=FALSE, cache=F----------------------
#  LDL_Vrnc_Chng_HDL <- LDL_HDLVrnc(dfSmplAHDLChngCV,
#                                              dfSmplA$CHOL,
#                                              dfSmplA$TG,
#                                              bootStrpReps = 2000)

## ---- fig.height=3, fig.width=8, position="!h", eval=FALSE,fig.cap="Variance of LDL (error propagation and boostrap for 2000 iterations), when the HDL variance changes."----
#  LDL_Vrnc_Chng_HDL <- do.call(cbind.data.frame, LDL_Vrnc_Chng_HDL)
#  
#  # Create a data frame with the changing HDL variance and the error propagation
#  # and bootstrap variances as columns, in order to use it for ggplot and
#  # remove the original data frame.
#  dfHDLChngVrnc <- cbind.data.frame(HDLVrnc= apply(dfSmplAHDLChngCV, 2, var),
#                   LDLErrPrp = LDL_Vrnc_Chng_HDL$ErrPropVrnc,
#                   LDLBootstrp = LDL_Vrnc_Chng_HDL$BootVrnc)
#  rm(LDL_Vrnc_Chng_HDL)
#  
#  dfHDLChngVrnc <- pivot_longer(dfHDLChngVrnc, cols=2:3,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  pltLDL_HDLVrnc <-
#    ggplot(data=dfHDLChngVrnc,
#           aes(x=HDLVrnc, y=Variance, linetype=Variance_Type)) +
#      geom_line() +
#      xlab("HDL Variance") + ylab("LDL Variance") +
#      scale_linetype_discrete(name="Variance Type",
#                              labels = c("Bootstrap Variance",
#                                         "Error Propagation Variance"))
#  pltLDL_HDLVrnc

## -----------------------------------------------------------------------------
CV(dfSmplA$TG)

## ---- TGVrncChng, fig.height=6, fig.width=8, fig.align="center", fig.cap="Means and variances for the TG distributions, in order to check that the mean stays constant and the variance increases.", cache=F----
dfSmplATGChngCV <- CV_Range(dfSmplA$TG,
                                       lower_CV_Bound = 0,
                                       upper_CV_Bound = 20,
                                       maxRandIter = 10000,
                                       plot=T)

## ---- LDL_TGVrnc2000,eval=FALSE,results='hide',cache=F------------------------
#  LDL_Vrnc_Chng_TG <- LDL_TGVrnc(dfSmplATGChngCV,
#                                              dfSmplA$CHOL,
#                                              dfSmplA$HDL,
#                                              bootStrpReps = 2000)

## ---- fig.height=3,eval=FALSE, fig.width=8, fig.cap="Variance of LDL (error propagation and boostrap for 2000 iterations), when the TG variance changes."----
#  LDL_Vrnc_Chng_TG <- do.call(cbind.data.frame, LDL_Vrnc_Chng_TG)
#  
#  # Create a data frame with the changing TG varince and the error propagation
#  # and bootstrap variances as columns, in order to use it for ggplot and
#  # remove the original data frame.
#  dfTGChngVrnc <- cbind.data.frame(TGVrnc= apply(dfSmplATGChngCV, 2, var),
#                   LDLErrPrp = LDL_Vrnc_Chng_TG$ErrPropVrnc,
#                   LDLBootstrp = LDL_Vrnc_Chng_TG$BootVrnc)
#  rm(LDL_Vrnc_Chng_TG)
#  
#  dfTGChngVrnc <- pivot_longer(dfTGChngVrnc, cols=2:3,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  pltLDL_TGVrnc <-
#    ggplot(data=dfTGChngVrnc,
#           aes(x=TGVrnc, y=Variance, linetype=Variance_Type)) +
#      geom_line() +
#      xlab("TG Variance") + ylab("LDL Variance") +
#      scale_linetype_discrete(name="Variance Type",
#                              labels = c("Bootstrap Variance",
#                                         "Error Propagation Variance"))
#  pltLDL_TGVrnc

## -----------------------------------------------------------------------------
CV(dfSmplA$HDL)

## ---- HDLVrncChng_AIP, fig.height=6, fig.width=8, fig.align="center", fig.cap="Means and variances for the HDL distributions, in order to check that the mean stays constant and the variance increases.", cache=F----
dfSmplAHDLChngCV <- CV_Range(dfSmplA$HDL,
                                       lower_CV_Bound = 0,
                                       upper_CV_Bound = 20,
                                       maxRandIter = 10000,
                                       plot=T)

## ---- AIP_HDLVrnc2000, eval=FALSE,cache=F, results='hide'---------------------
#  AIP_Vrnc_Chng_HDL <- AIP_HDLVrnc(dfSmplAHDLChngCV,
#                                            dfSmplA$TG,
#                                            SI=FALSE,
#                                            bootStrpReps = 2000)

## ---- fig.height=3,  eval=FALSE, fig.width=8, fig.cap="Variance of AIP (error propagation and boostrap for 2000 iterations), when the HDL variance changes."----
#  AIP_Vrnc_Chng_HDL <- do.call(cbind.data.frame, AIP_Vrnc_Chng_HDL)
#  
#  # Create a data frame with the changing HDL varince and the error propagation
#  # and bootstrap variances as columns, in order to use it for ggplot and
#  # remove the original data frame.
#  dfHDLChngVrnc <- cbind.data.frame(HDLVrnc= apply(dfSmplAHDLChngCV, 2, var),
#                   AIPErrPrp = AIP_Vrnc_Chng_HDL$ErrPropVrnc,
#                   AIPErrPrp2Ord = AIP_Vrnc_Chng_HDL$ErrPropVrnc2Ord,
#                   AIPBootstrp = AIP_Vrnc_Chng_HDL$BootVrnc)
#  rm(AIP_Vrnc_Chng_HDL)
#  
#  dfHDLChngVrnc <- pivot_longer(dfHDLChngVrnc, cols=2:4,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  pltAIP_HDLVrnc <-
#    ggplot(data=dfHDLChngVrnc,
#           aes(x=HDLVrnc, y=Variance, linetype=Variance_Type)) +
#      geom_line() +
#      xlab("HDL Variance") + ylab("AIP Variance") +
#      scale_linetype_discrete(name="Variance Type",
#                              labels = c("Bootstrap Variance",
#                                         "Error Propagation Variance",
#                                         "Error Propagation Variance 2nd Order"))
#  pltAIP_HDLVrnc

## -----------------------------------------------------------------------------
CV(dfSmplA$TG)

## ---- TGVrncChng_AIP, fig.height=6, fig.width=8, fig.align="center", fig.cap="Means and variances for the TG distributions, in order to check that the mean stays constant and the variance increases.", cache=F----
dfSmplATGChngCV <- CV_Range(dfSmplA$TG,
                                       lower_CV_Bound = 0,
                                       upper_CV_Bound = 20,
                                       maxRandIter = 10000,
                                       plot=T)

## ---- AIP_TGVrnc2000,eval=FALSE,results='hide', cache=F-----------------------
#  AIP_Vrnc_Chng_TG <- AIP_TGVrnc(dfSmplATGChngCV,
#                                            dfSmplA$HDL,
#                                            SI=FALSE,
#                                            bootStrpReps = 2000)

## ---- fig.height=3,eval=FALSE, fig.width=8, fig.cap="Variance of AIP (error propagation and boostrap for 2000 iterations), when the TG variance changes."----
#  AIP_Vrnc_Chng_TG <- do.call(cbind.data.frame, AIP_Vrnc_Chng_TG)
#  
#  # Create a data frame with the changing TG variance and the error propagation
#  # and bootstrap variances as columns, in order to use it for ggplot and
#  # remove the original data frame.
#  dfTGChngVrnc <- cbind.data.frame(TGVrnc= apply(dfSmplATGChngCV, 2, var),
#                   AIPErrPrp = AIP_Vrnc_Chng_TG$ErrPropVrnc,
#                   AIPErrPrp2Ord = AIP_Vrnc_Chng_TG$ErrPropVrnc2Ord,
#                   AIPBootstrp = AIP_Vrnc_Chng_TG$BootVrnc)
#  rm(AIP_Vrnc_Chng_TG)
#  
#  dfTGChngVrnc <- pivot_longer(dfTGChngVrnc, cols=2:4,
#                               names_to = "Variance_Type",
#                               values_to = "Variance")
#  pltAIP_TGVrnc <-
#    ggplot(data=dfTGChngVrnc,
#           aes(x=TGVrnc, y=Variance, linetype=Variance_Type)) +
#      geom_line() +
#      xlab("TG Variance") + ylab("AIP Variance") +
#      scale_linetype_discrete(name="Variance Type",
#                              labels = c("Bootstrap Variance",
#                                         "Error Propagation Variance",
#                                         "Error Propagation Variance 2nd Order"))
#  pltAIP_TGVrnc

## ----ggplot2, warning = FALSE,eval=FALSE,echo=FALSE, out.width="50%"----------
#  pltLDL_CHOLVrnc
#  pltLDL_HDLVrnc
#  pltLDL_TGVrnc
#  pltAIP_HDLVrnc
#  pltAIP_TGVrnc

