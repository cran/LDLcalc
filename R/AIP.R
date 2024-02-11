#' @title Calculation of Atherogenic Index of Plasma (AIP)
#' @description Calculates the Atherogenic Index of Plasma (AIP) from triglyceride and HDL values.
#' @param TG A vector or data frame column containing the triglyceride (TG) values to be used for the
#' calculation of the atherogenic index of plasma (AIP). TG and HDL must be of the same length.
#' @param HDL A vector or data frame column containing the high density lipoprotein (HDL) values
#' to be used for the calculation of the atherogenic index of plasma (AIP). TG and HDL
#' must be of the same length.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI units for TG and HDL
#' (mmol/L). If mg/dl units are provided instead, SI must be set to FALSE.
#' @param roundDigit Decimal digits to round the result to (default = 4).
#' @usage AIPcalc(TG, HDL, SI ,roundDigit)
#' @return A vector of AIP values of length equal to the length of TG and HDL.
#' @references Dobiasova M, Frohlich J. The plasma parameter log (TG/HDL-C)
#' as an atherogenic index: correlation with lipoprotein particle size and
#' esterification rate in apo B- lipoprotein-depleted plasma (FERHDL).
#' Clin Biochem. 2001, 34:583-88.
#'
#' Tan MH, Johns D, Glazer NB. Pioglitazone reduces atherogenic index of plasma
#' in patients with type 2 diabetes. Clin Chem. 2004, 50:1184-88.
#'
#'  Nwagha UI, Ingweh JC. A significant indicator for the onset of
#' atherosclerosis during menopause in hypertensive females of South East Nigeria.
#' J Coll Med. 2005, 10(2):67-71.
#'
#' Daniels LB, Laughlin G, Sarno MJ. Lp- PLA2 is an independent predictor of
#' incidence of coronary heart disease in apparently healthy older population.
#' J Am Col Cardiol. 2008,51:913-19.
#'
#' Xiaowei Z, Lugang Y, Hui Z, Qinhua M, Xiaohua Z, Ting L, et al.
#' Atherogenic index of plasma is a novel and better biomarker associated with
#' obesity: A population- based cross-sectional study in China.
#' Lipids Health Dis. 2018,17(1):37.
#' @examples
#' \dontrun{
#' AIP = AIPcalc(sampleA$TG, sampleA$HDL)
#' }
#' @export
AIPcalc<-function(TG, HDL, SI=TRUE, roundDigit=4){
  if ( length(TG)!=length(HDL)){
    message("The length of TG and HDL the must be the same")
    stopifnot()}
    if (SI==FALSE){
      TG=TG*0.01129
      HDL=HDL*0.0259
    }
      AIPcalcul <- round(log10((TG)/(HDL)), roundDigit)
    return(AIPcalcul)
  }


#' @title Variance of Atherogenic Index of Plasma (AIP) using Error Propagation.
#' @description Calculate the variance of AIP using Error Propagation (The Delta Method)
#' @param TG A vector or data frame column containing the triglyceride (TG) values to be used for
#' the calculation of the variance of Atherogenic Index of Plasma (AIP).
#' TG and HDL must be of the same length.
#' @param HDL A vector or data frame column containing the high density lipoprotein (HDL) values
#' to be used for the calculation of variance of the Atherogenic Index of Plasma (AIP).
#' TG and HDL must be of the same length.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI units for TG and HDL
#' (mmol/L). If mg/dl units are provided instead, SI must be set to FALSE.
#' @param roundDigit Decimal digits to round the result to (default = 5).
#' @return The variance of AIP using error propagation theory.
#' @references Casella G, Berger RL. Statistical Inference. 2nd ed.
#' Duxbury Thomson Learning, 2002, Pages 240-245
#'
#' Joint Commitee for Guides in Metrology. Evaluation of measurement data-
#' Guide to the expression of uncertainty in measurement. 2008.
#'
#' Joint Commitee for Guides in Metrology. International Vocabulary of Metrology (VIM)-
#' Basic and General Concepts and Associated Terms. 2012.
#' @examples
#' \dontrun{
#' AIPpropagationVar = AIPErrPrp(sampleA$TG,sampleA$HDL)
#' }
#' @export
AIPErrPrp=function(TG, HDL, SI=TRUE, roundDigit=5) {
  if (SI==FALSE){
    TG=TG*0.01129
    HDL=HDL*0.0259
  }
  AIPpropVar=round(fAIPErrPrp(TG, HDL),roundDigit)
  return(AIPpropVar)
}


#' @title Second Order Taylor expansion AIP Error propagation variance
#' @description Calculate the second order Taylor Expansion error propagation
#' (delta method) variance of the Atherogenic Index of Plasma.
#' @param TG A vector or data frame column containing the triglyceride (TG) values to be used for
#' the calculation of the variance of the Atherogenic Index of Plasma (AIP)
#' using second order Taylor expansion. TG and HDL must be of the same length.
#' @param HDL A vector or data frame column containing the high density lipoprotein (HDL) values
#' to be used for the calculation of the variance of the Atherogenic Index of Plasma (AIP)
#' using second order Taylor expansion.
#' TG and HDL must be of the same length.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI units for TG and HDL
#' (mmol/L). If mg/dl units are provided instead, SI must be set to FALSE.
#' @param roundDigit Decimal digits to round the result to (default = 5).
#' @usage AIPErrPrp2Ord(TG, HDL, SI=TRUE, roundDigit=5)
#' @return The variance of AIP using error propagation theory.
#' @references Casella G, Berger RL. Statistical Inference. 2nd ed.
#' Duxbury Thomson Learning 2002, Pages 240-245
#'
#' Joint Commitee for Guides in Metrology. Evaluation of measurement data-
#' Guide to the expression of uncertainty in measurement. 2008.
#'
#' Joint Commitee for Guides in Metrology. International Vocabulary of Metrology (VIM)-
#' Basic and General Concepts and Associated Terms. 2012.
#'
#' @examples
#' \dontrun{
#' AIPpropagationVarTaylor = AIPErrPrp2Ord(sampleA$TG,sampleA$HDL)
#' }
#' @export
AIPErrPrp2Ord <-function(TG, HDL, SI=TRUE, roundDigit=5) {
  if (SI==FALSE){
    TG=TG*0.01129
    HDL=HDL*0.0259
  }
  AIPpropVarTaylor=round(fAIPErrPrp2Ord(TG, HDL), roundDigit)
  return(AIPpropVarTaylor)
  }


#' @title Calculate variance of AIP using Bootstrapping.
#' @description  This function can be used to calculate the variance of the
#' Atherogenic Index of Plasma (AIP) using Bootstrapping.
#' @param TG A vector or data frame column containing the triglyceride (TG)
#' values to be used for
#' the calculation of the variance of the Atherogenic Index of Plasma (AIP)
#' using bootstrapping. TG and HDL must be of the same length.
#' @param HDL A vector or data frame column containing the high density
#' lipoprotein (HDL) values
#' to be used for the calculation of the variance of the Atherogenic Index
#' of Plasma (AIP)
#' using bootstrapping.
#' TG and HDL must be of the same length.
#' @param sampleSize (default = length of TG or HDL) The sample size that will be
#' generated at each bootstrapping sampling
#' round. Size of bootstrapped samples cannot be larger than the original.
#' @param SI Boolean (default=TRUE). AIP is by definition calculated using SI
#' units for TG and HDL
#' (mmol/L). If mg/dl units are provided instead, SI must be set to FALSE.
#' @param noOfReps  (default = 1000) Number of repetitions of the bootstrapping.
#' @param pb Display a progress bar (default = FALSE)
#' @return It returns a data table with four columns. The first
#' column contains the mean of the AIP values for each iteration. The second
#' column contains the median of each iteration.The third
#' column contains the variance and the fourth column contains
#' the CV of each iteration. It also returns the median of the "Mean", "Var"
#' and "CV" columns of the data table.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' AIPbootstrVar = AIPbootVrnc(sampleA$TG,sampleA$HDL)
#' }
#' @export
AIPbootVrnc <- function(TG, HDL, sampleSize=length(TG), SI=TRUE, noOfReps=1000, pb=F) {
  # Μετατροπή σε SI μονάδες αν απαιτείται
  if (SI == FALSE) {
    TG <- TG * 0.01129
    HDL <- HDL * 0.0259
  }
  
  # Έλεγχος για ίσο μήκος παραμέτρων
  if((length(HDL) != length(TG))) {
    message("The two parameters must have the same number of measurements")
    return(NULL)
  }
  
  n <- length(TG)
  if(sampleSize > n) {
    message("Size of bootstrapped datasets cannot be larger than the original dataset. Setting sampleSize to n")
    sampleSize <- n
  }
  
  dfTmp <- data.frame(TG, HDL)
  
  # Αρχικοποίηση λίστας για αποθήκευση των αποτελεσμάτων
  resultsList <- list()
  
  if(pb) pb <- txtProgressBar(min = 0, max = noOfReps, style = 3)
  
  for(bootIdx in 1:noOfReps) {
    dfSmpl <- dfTmp[sample(nrow(dfTmp), size=sampleSize, replace = TRUE), ]
    dfAIPSmpl <- log10(dfSmpl$TG / dfSmpl$HDL)
    AIPSmplMean <- mean(dfAIPSmpl)
    AIPSmplMedian <- median(dfAIPSmpl)
    AIPSmplVar <- var(dfAIPSmpl)
    AIPSmplCV <- CV(dfAIPSmpl) # Υποθέτουμε ότι η συνάρτηση CV είναι ήδη ορισμένη
    
    # Προσθήκη των αποτελεσμάτων στη λίστα
    resultsList[[bootIdx]] <- list(Mean = AIPSmplMean, Median = AIPSmplMedian, Var = AIPSmplVar, CV = AIPSmplCV)
    
    if(pb) setTxtProgressBar(pb, bootIdx)
  }
  
  # Μετατροπή της λίστας αποτελεσμάτων σε data.table
  dtAIPBoot <- data.table::rbindlist(resultsList)
  
  # Ορισμός ονομάτων στηλών
  stats::setNames(dtAIPBoot, c("Mean", "Median", "Var", "CV"))
  
  return(list("dataTable" = dtAIPBoot, "Mean" = median(dtAIPBoot$Mean),
              "Var" = median(dtAIPBoot$Var), "CV" = median(dtAIPBoot$CV)))
}