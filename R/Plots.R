#' @title Plot variance density
#' @description Plot density of the Variance along with the median, 2.5 and 97.5
#' percentile as vertical lines
#' @param Vector A vector or data frame columns with the values whose
#' @param title Title of the plot (default="").
#' distribution density is to be plotted.
#' @return The function returns no values but plots the density plot.
#' @examples
#' \dontrun{
#' LDLbootstrVar=as.data.frame(LDLbootVrnc(sampleA$CHOL,sampleA$HDL,sampleA$TG))
#' DensityPlotOfVar(LDLbootstrVar$dataTable.CV)
#' }
#' @export
DensityPlotOfVar <- function(Vector, title="") {
  plot(density(Vector), main=title)
  abline(v=median(Vector))
  v1=quantile(Vector, probs=c(0.025))
  abline(v=v1, lty=2)
  v2=quantile(Vector, probs=c(0.975))
  abline(v=v2, lty=2)
 return(list(v1,v2))
}

#' @title Plot discrete histogram
#' @description Plot discrete histogram (barplot) of the data frame column named in the param
#' argument. It also plots the mean as a vertical continuous line, the mean plus/minus 2
#' standard deviations as veritcal dotted lines and overlays a density plot of the
#' normal distribution with mean and standard deviation corresponding to those of the data.
#' @param DF A data frame containing columns whose discrete histogram is to be plotted.
#' @param param The column name of the columns for which to plot the discrete histogram.
#' @param title Title of the plot (default="").
#' @importFrom stats sd dnorm
#' @return It returns a ggplot object.
#' @examples
#' \dontrun{
#' PlotDiscrHist(sampleA,"LDL")
#' }
#' @export
PlotDiscrHist <- function(DF, param, title="") {
  df = data.frame(DF)
  plt <-
    ggplot2::ggplot(data=df, ggplot2::aes_string(x=param)) + ggplot2::geom_bar() +
    ggplot2::stat_function(fun=fScaleNormDistr, args=c(mean=mean(df[,param]), sd=sd(df[,param]), n=nrow(df))) +
    ggplot2::geom_vline(xintercept = mean(df[,param]), size=1) +
    ggplot2::geom_vline(xintercept = mean(df[,param]) - 2*sd(df[,param]), linetype=2) +
    ggplot2::geom_vline(xintercept = mean(df[,param]) + 2*sd(df[,param]), linetype=2) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(plt)
}

#' @title Scatterplot of pairs of parameters
#' @description Scatterplot of pairs of parameters with the corresponding regression line.
#' @param df The data frame with the parameterrs to be plotted
#' @param xParam The parameter (column name) to be plotted in the abscissa (x axis).
#' @param yParam The prameter (column name) to be plotted in the ordinate (y axis).
#' @return The function returns a ggplot2 object.
#' @importFrom ggplot2 ggplot geom_point geom_smooth theme ggtitle
#' @examples
#' \dontrun{
#' PlotCorrWithRegrLine(sampleA,"CHOL", "HDL")
#' }
#' @export
  PlotCorrWithRegrLine <- function(df,xParam, yParam) {
  plt <-
    ggplot(data=df, ggplot2::aes_string(x=xParam, y=yParam)) +
    geom_point() +
    geom_smooth(method="lm", se=F) +
    ggtitle(paste(xParam, " versus ", yParam, sep="")) +
    theme(plot.title = ggplot2::element_text(hjust = 0.5))
  return(plt)
}


#' @title Plot density of bootstrap and empirical variance for LDL
#' @description  Plot density of the bootstrap variance for LDL along with the median,
#' 2.5 and 97.5 percentile as vertical lines.
#' Also plots the empirical variance as a segment from the bottom to the middle
#' of the graph and the error propagation variance as a line segment from the
#' top to the middle of the graph. DF must the dataframe of the list
#' returned from the respective Bootstrap variance function.
#' @param DF The data frame containing the bootstrap variances of LDL.
#' @param title Title of the plot (default="")
#' @param empirVrnc The value of of the empirical (experimental) variance.
#' @param errPropVrnc The value of the error propagation variance.
#' @return Plots the respective plot.
#' @examples
#' \dontrun{
#' LDL_empirVrnc = var(sampleA$LDL)
#' LDL_errPropVrnc = LDLErrPrp(sampleA$CHOL,sampleA$HDL,sampleA$TG)
#' LDLbootStrp=as.data.frame(LDLbootVrnc(sampleA$CHOL,sampleA$HDL,sampleA$TG))
#' LDL_DensityPlotOfbootst(LDLbootStrp,"Title",LDL_empirVrnc,LDL_errPropVrnc)
#' }
#' @export
LDL_DensityPlotOfbootst <- function(DF, title="", empirVrnc, errPropVrnc) {
  xMin <- min(density(DF$Var)$x)
  xMax <- max(density(DF$Var)$x)
  yMin <- min(density(DF$Var)$y)
  yMax <- max(density(DF$Var)$y)
  plot(density(DF$Var), main=title, xlab="LDL Variance", xlim=c(0,xMax))
  abline(v=median(DF$Var))
  abline(v=quantile(DF$Var, probs=c(0.025)), lty=2)
  abline(v=quantile(DF$Var, probs=c(0.975)), lty=2)
  segments(x0=empirVrnc,x1=empirVrnc, y0=yMin, y1=yMax/2, lty=3, lwd=2)# Plot segment corresponding to
  #empirical variance
  segments(x0=errPropVrnc, x1=errPropVrnc, y0=yMax, y1=yMax/2, lty=4, lwd=2)# Plot segment corresponding
  # to error propagation variance
}

#' @title Plot density of bootstrap and empirical variance for AIP
#' @description  Plot density of the bootstrap variance for AIP along with the median,
#' 2.5 and 97.5 percentile as vertical lines. Also plots the empirical variance
#' as a segment from the bottom to the middle of the graph and the error
#' propagation variance as a line segment from the top to the middle of the
#' graph. DF must the dataframe of the list returned from the respective
#' (LDL or AIP) Bootstrap variance function.
#' @param DF A data frame containing the bootstrap variances of AIP.
#' @param title Title of the plot (default="").
#' @param empirVrnc Value of the empirical (experimental) variance of AIP.
#' @param errPropVrnc Value of the first order error propagation variance of AIP.
#' @param errPropVrnc2Ord Value of the second order error propagation variance of AIP.
#' @return Plots the respective plot.
#' @importFrom graphics par abline segments
#' @importFrom stats density median quantile
#' @examples
#' \dontrun{
#' sampleA$AIP = AIPcalc(sampleA$TG,sampleA$HDL, SI=FALSE)
#' AIP_empirVrnc=var(sampleA$AIP)
#' AIP_errPropVrnc=AIPErrPrp(sampleA$TG,sampleA$HDL, SI=FALSE)
#' AIP_errPropVrnc2Ord=AIPErrPrp2Ord(sampleA$TG,sampleA$HDL, SI=FALSE)
#' DfAIPboost=as.data.frame(AIPbootVrnc(sampleA$TG,sampleA$HDL, SI=FALSE))
#' AIP_DensityPlotOfbootst(DfAIPboost,"Title",AIP_empirVrnc, AIP_errPropVrnc, AIP_errPropVrnc2Ord)
#' }
#' @export
AIP_DensityPlotOfbootst <- function(DF, title="", empirVrnc, errPropVrnc, errPropVrnc2Ord) {
  xMin <- min(density(DF$Var)$x)
  xMax <- max(density(DF$Var)$x)
  yMin <- min(density(DF$Var)$y)
  yMax <- max(density(DF$Var)$y)
  plot(density(DF$Var), main=title, xlim=c(0,xMax), xlab="AIP Variance")
  abline(v=median(DF$Var))
  abline(v=quantile(DF$Var, probs=c(0.025)), lty=2)
  abline(v=quantile(DF$Var, probs=c(0.975)), lty=2)
  segments(x0=empirVrnc,x1=empirVrnc, y0=yMin, y1=yMax/2, lty=3, lwd=2)# Plot segment corresponding
  # to empirical variance
  segments(x0=errPropVrnc, x1=errPropVrnc, y0=yMax *0.25, y1=yMax*0.75, lty=5, lwd=2)# Plot segment corresponding
  # to error propagation variance
  segments(x0=errPropVrnc2Ord, x1=errPropVrnc2Ord, y0=yMax, y1=yMax/2, lty=4, lwd=2)# Plot segment
  # corresponding to 2nd Order error propagation variance
}



