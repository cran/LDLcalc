##flst2DF - Convert a list (e.g. from lapply) to a data frame.
flst2DF <- function(lst) {
  DF = as.data.frame(do.call(rbind, lst))
}

`%tin%` <- function(x, y) {
  mapply(assign, as.character(substitute(x)[-1]), y,
         MoreArgs = list(envir = parent.frame()))
  invisible()
}

#fScaleNormDistr
#Scale normal distribution frequencies by a factor of n, so as to add the normal
# curve to a barplot (discrete histogram).
#It is to be used in fPlotDiscrHist
fScaleNormDistr <- function(x, mean, sd, n) {
  dnorm(x = x, mean = mean, sd = sd) * n
}

## fAIPErrPrp
fAIPErrPrp <- function(TG, HDL) {
  (  mean(TG*0.0113)^2 * var(HDL*0.0259)   +
       mean(HDL*0.0259)^2 * var(TG*0.0113)   -
       2 *  mean(TG*0.0113) * mean(HDL*0.0259) * cov(TG*0.0113,HDL*0.0259)
  )   /
    (log(10)^2 * mean(TG*0.0113)^2 * mean(HDL*0.0259)^2)
}


# fAIPErrPrp2Ord
# x1=TG, x2=HDL
fAIPErrPrp2Ord <- function(x1,x2) {
  u1 <- mean(x1)
  u2 <- mean(x2)
  m2_1 <- moments::moment(x1,order=2,central=T)#4th central moment of x1
  m2_2 <- moments::moment(x2,order=2,central=T)#4th central moment of x2
  m4_1 <- moments::moment(x1,order=4,central=T)#4th central moment of x1
  m4_2 <- moments::moment(x2,order=4,central=T)#4th central moment of x2
  corx1x2 <- cor(x1, x2)
  nominator <-
    u2^4 * m4_1 -
    4 * u1^3 * u2^3 * corx1x2 * sd(x1) * sd(x2) +
    2 * u1^2 * u2^2 * m2_1 * (u2 - corx1x2 * sd(x2)) * (u2 + corx1x2 * sd(x2)) +
    u1^4 * (2 * u2^2 * m2_2 + m4_2)

  denom <- 2 * u1^4 * u2^4 * log(10)^2

  AIPVrnc2Ord <- nominator / denom
  return(AIPVrnc2Ord)
}

#fScaleNormDistr
#Scale normal distribution frequencies by a factor of n, so as to add the normal curve to a barplot (discrete histogram).
# It is to be used in fPlotDiscrHist

fScaleNormDistr <- function(x, mean, sd, n) {
  + dnorm(x = x, mean = mean, sd = sd) * n
}

#exclude_runif
#Function to give a random integer number using the random uniform distribution after excluding a specified value.
#' @importFrom stats runif
exclude_runif <- function(n, min, max, valueToExclude) {
  res = valueToExclude
  while (res == valueToExclude) {
    res = round(runif(1, min, max))
  }
  return(res)
}

fPltDensMedPerc <- function(DF, title) {
  plot(density(DF$dataTable$Var), main=title)
  abline(v=median(DF$dataTable$Var))
  abline(v=quantile(DF$dataTable$Var, probs=c(0.025)), lty=2)
  abline(v=quantile(DF$dataTable$Var, probs=c(0.975)), lty=2)
}

#fChngDistVrncMeanConst
#Function to take a set of values and increase its variance while keeping the sample mean constant
fChngDistVrncMeanConst <- function(grp) {
  grpLen = length(grp)
  indxOfValueToDecrease = round(runif(1, 1, grpLen))
  indxOfValueToIncrease = exclude_runif(1, 0, grpLen, indxOfValueToDecrease )
  grpNew = grp
  grpNew[indxOfValueToDecrease] = grpNew[indxOfValueToDecrease] - 1
  grpNew[indxOfValueToIncrease] = grpNew[indxOfValueToIncrease] + 1
  return(grpNew)
}

#fIncrCVMeanConst
#Function to take a set of values and increase its CV (using function fChngDistVrncMeanConst) until the
# CV reaches the maximum specified value.
fIncrCVMeanConst <- function(grp, maxCV, maxIter=10000) {
  grpInit = vector(mode="numeric", length = length(grp))
  grpCombIncrVariance = data.frame(cbind(grpInit, grp))
  currCV = CV(grp)
  iterCounter = 0
  while( (currCV<=maxCV & iterCounter<maxIter)) {
    grpOld = grpCombIncrVariance[,ncol(grpCombIncrVariance)]
    grpNew = fChngDistVrncMeanConst(grpOld)
    if( (var(grpNew)>var(grpOld)) & (mean(grpNew)==mean(grpOld)) ) {
      grpCombIncrVariance = cbind( grpCombIncrVariance,grpNew)
    }
    currCV = CV(grpCombIncrVariance[,ncol(grpCombIncrVariance)])

  }
  grpCombIncrVariance$grpInit = NULL
  return(grpCombIncrVariance)

}


##fDecrCVMeanConst
#Function to take a set of values and decrease its CV (using function fChngDistVrncMeanConst) until the CV
#reaches the minimun specified value.
fDecrCVMeanConst <- function(grp, minCV, maxIter=10000) {
  iterCounter = 1
  grpInit = vector(mode="numeric", length = length(grp))
  grpCombDecVariance = data.frame(cbind(grpInit, grp))
  currCV = CV(grp)
  while( (currCV>minCV) & (iterCounter<maxIter) ) {
    grpOld = grpCombDecVariance[,ncol(grpCombDecVariance)]
    grpNew = fChngDistVrncMeanConst(grpOld)
    if( (var(grpNew)<var(grpOld)) & (mean(grpNew)==mean(grpOld)) ) {
      grpCombDecVariance = cbind( grpCombDecVariance, grpNew)
    }
    currCV = CV(grpCombDecVariance[,ncol(grpCombDecVariance)])
    iterCounter = iterCounter + 1
  }
  grpCombDecVariance$grpInit = NULL
  return(grpCombDecVariance)
}

#fCombineVariances
#Function to combine the two data frames that resulted from fIncrCVMeanConst() and fDecrCVMeanConst() into one data frame with the samples arranged in increasing variance order.
fCombineVariances <- function(dfVariancesDec, dfVariancesIncr, plot=F) {
  dfVariancesDecRev = rev(dfVariancesDec)

  dfVariances = cbind(dfVariancesDecRev, dfVariancesIncr[,2:ncol(dfVariancesIncr)])
  if(plot) {
    par(mfrow=c(2,1))
    plot(1:ncol(dfVariances), apply(dfVariances, 2, mean), type="l",xlab = "Modified Samples", ylab = "Sample Mean",
         main ="Results of CVIncrDecr")
    plot(1:ncol(dfVariances), apply(dfVariances, 2, CV), xlab = "Modified Samples", ylab = "Sample CV")
    par(mfrow=c(1,1))
  }
  return(dfVariances)
}

#fChngDistrCVMeanConst
# Function that combines most of the previous functions in order to change the CV of a distribution while
#keeping the mean constant. An upper and a lower value is specified for the CV.
# It uses the functions fIncrCVMeanConst(), fDecrCVMeanConst() and fCombineVariances(), returning a DF with the columns of samples
# of different CVs in ascending order and also plots the change of mean and cv values of the sample.
fChngDistrCVMeanConst <- function(vecDistr, lowerCV, upperCV, maxIter=10000, plot=F) {
  dfDecrCV = fDecrCVMeanConst(vecDistr, minCV=lowerCV, maxIter = maxIter)
  dfIncrCV = fIncrCVMeanConst(vecDistr, maxCV=upperCV, maxIter = maxIter)
  dfCombVrncs = fCombineVariances(dfDecrCV, dfIncrCV, plot=plot)
  return(dfCombVrncs)
}

fAverageData <- function(dat, sampleSize, noOfReps) {
  vecMeans = vector(mode="numeric", length=noOfReps)
  10
  for(i in seq(1,noOfReps)) {
    smpl = sample(dat, size=sampleSize, replace = T)
    vecMeans[i] = mean(smpl)
  }
  return(vecMeans)
}
