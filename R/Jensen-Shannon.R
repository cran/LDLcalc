# This .R file requires the "janitor" and "philentropy" packages to run.

#' @title Calculate the Jensen-Shannon divergence (JSD)
#' @description Calculate the Jensen-Shannon divergence between the same
#' parameter of different
#' data sets. This works only for two distributions.
#' @param vec1 The vector containing the values of the first data set.
#' @param vec2 The vector containing the values of the second data set.
#' @return The function returns the value of the Jensen-Shannon divergence (JSD),
#' the frequencies of the two datasets (freq1, freq2, dfFreqs) and the
#' probabilities of the frequencies of the two groups, which should sum up
#' to 1.
#' @references D.M. Endres, J.E. Schindelin, A new metric for probability distributions, IEEE
#' Trans. Inf. Theory (2003), https://doi.org/10.1109/TIT.2003.813506.
#'
#' F. Oesterreicher, I. Vajda, A new class of metric divergences on probability spaces
#' and its applicability in statistics, Ann. Inst. Stat. Math. (2003), https://doi.org/
#'  10.1007/BF02517812.
#' @examples
#' \dontrun{
#' JSDCalc = JSD(model[["trainData"]]$LDLd,model[["testData"]]$LDLd)
#' }
#' @export
  JSD <- function(vec1, vec2) {
  #Calculate the frequencies of vec1
  df1Freqs <- as.data.frame(janitor::tabyl(vec1))
  #The tabyl command gives a data frame with not only the frequencies but also the
  #percentages. Remove this column
  df1Freqs$percent <- NULL
  #Change column names
  colnames(df1Freqs) <- c("Param", "vec1Freq")
  #The same as above for vec2
  df2Freqs <- as.data.frame(janitor::tabyl(vec2))
  df2Freqs$percent <- NULL
  colnames(df2Freqs) <- c("Param", "vec2Freq")

  #Merge the two frequency data frames by the "Param" column.
  #When a Param value does not exist in one of the data frames, an NA is
  #inserted
  dfFreqs <- merge(df1Freqs, df2Freqs, by="Param", all.x=T, all.y = T)
  colnames(dfFreqs) <- c("Param", "vec1", "vec2")
  dfFreqs[is.na(dfFreqs)] <- 0 #Replace NA with zeros.

  #Calc the probabilities of vec1 and vec2 frequencies. They should sum to 1.
  dfFreqs$vec1Prob <- dfFreqs$vec1 / sum(dfFreqs$vec1)
  vec1ProbSum <- sum(dfFreqs$vec1Prob)
  dfFreqs$vec2Prob <- dfFreqs$vec2 / sum(dfFreqs$vec2)
  vec2ProbSum <- sum(dfFreqs$vec2Prob)

  JSDres <-philentropy::JSD(rbind(dfFreqs$vec1Prob, dfFreqs$vec2Prob))
  return(list(JSD=JSDres, freq1=df1Freqs, freq2=df2Freqs,
              freqs=dfFreqs, vec1PSum=vec1ProbSum, vec2PSum=vec2ProbSum))
}


#' @title Calculate the Jensen-Shannon divergence (JSD) between a discrete
#' empirical distribution and the normal distribution.
#' @description Calculates the Jensen-Shannon divergence between a discrete
#' distribution and the corresponding normal distribution with mean and standard
#' deviation the same as these of the discrete one.
#' @param dfSmpl A data frame containing the values of the discrete distribution.
#' The data frame may contain more that one column with discrete distribution
#' values. The argument "param" specified next will determine which column
#' will be used
#' @param param The name of the column to be used.
#' @return The function returns the Jensen-Shannon divergence between the
#' discrete and corresponding normal distribution. It also returns a data frame
#' with the empirical probability of the values supplied in the column
#' as well as the empirical probabilies one of the normal discrete distribution.
#' distribution.
#' @references D.M. Endres, J.E. Schindelin, A new metric for probability distributions, IEEE
#' Trans. Inf. Theory (2003), https://doi.org/10.1109/TIT.2003.813506.
#'
#' F. Oesterreicher, I. Vajda, A new class of metric divergences on probability spaces
#' and its applicability in statistics, Ann. Inst. Stat. Math. (2003), https://doi.org/
#'  10.1007/BF02517812.
#' @importFrom stats sd cov rnorm
#' @examples
#' \dontrun{
#' JSD.between.empirical.Normal =JSDNormal(sampleA,"LDL")
#' }
#' @export
JSDNormal <- function(dfSmpl, param) {
  dfSmplTbl = as.data.frame(table(dfSmpl[,param])) # Create a table of counts of the sample empirical distribution
  dfSmplTbl$prob = dfSmplTbl$Freq / sum(dfSmplTbl$Freq) # Calculate the probabilities
  d = round(rnorm(1e5, mean=mean(dfSmpl[,param]), sd=sd(dfSmpl[,param])), 0) #Create norm distr with mean,sd of the sample
  dfSmplNormTbl = as.data.frame(table(d)) # Calculate the frequencies (counts) of the created normal distr
  dfSmplNormTbl$prob = dfSmplNormTbl$Freq / sum(dfSmplNormTbl$Freq) #Calc probabilities of the created normal distr
  colnames(dfSmplTbl) = c("val", "freq", "prob") #Set colnames
  colnames(dfSmplNormTbl) = c("val", "freq", "prob")

  # Create df to hold the empirical and created (predicted) probabilities for each value of the empirical distr
  dfJSD = data.frame(probEmpirical=numeric(nrow(dfSmplTbl)), probPredicted=numeric(nrow(dfSmplTbl)))# Create dt
  counter=1
  #Loop to create a table with the probs of each value in the empirical distr and its corresponding prob in the created
  # normal distr
  for(idx in dfSmplTbl$val) { # For each value in the empirical distr
    if(idx %in% dfSmplNormTbl$val) { # if this value exists in the created distr
      dfJSD[counter, "probEmpirical"] = dfSmplTbl[counter,"prob"] # put in the df the empirical probability
      dfJSD[counter, "probPredicted"] = dfSmplNormTbl[which(dfSmplNormTbl$val==idx), c("prob")] # and the created one
    }
    else { # if the value in the empirical distr does not exist in the created normal distr
      dfJSD[counter,1] = dfSmplTbl[counter,"prob"] #put the empirical probability and the df
      dfJSD[counter,2] = 0 # and set the created (predicted) prpbability to 0.
    }
    counter = counter + 1
  }
  JSD = sqrt(philentropy::JSD(rbind(dfJSD$probEmpirical, dfJSD$probPredicted)))
  return(list("EmpirAndPredProb"=dfJSD, "JSD"=JSD))
}





#' @title Chebysev's inequality
#' @description Function to calculate lower and upper bounds for which at least 75% of data
#' points of a distribution lie within ±2SD of the mean using Chebyshev's inequality.

#' @param vec The vector containing the values of the  data set.
#' @return It outputs the lower and upper bound of Chebyshevs 75% range and the
#'  lower and upper value of the range of the observed (measured) distribution.
#'
#' @references Bienayme I. Considerations a l appui de la de couverte d
#'  laplace. Comptes Rendus de l Acade mie des Sciences 1853; 37: 309–324.
#'
#'. Chebyshev P. Des valeurs moyennes. Journal de Mathematiques Pures .
#'  Appliqee 1867; 2(12): 177–184s
#' @examples
#' \dontrun{
#' chebysevBounds = chebyshev(sampleA$LDL)
#' }
#' @export
chebyshev <- function(vec) {
  chebyshevLowerBound = round(mean(vec) - 2*sd(vec)) # Calculate lower bound
  chebyshevUpperBound = round(mean(vec) + 2*sd(vec)) # Calculate upper bound
  res = list(LowerBound= chebyshevLowerBound,
             UpperBound= chebyshevUpperBound,
             LowerRange= range(vec)[[1]],
             UpperRange = range(vec)[[2]])
  return(res)
}

