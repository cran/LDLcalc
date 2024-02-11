# This .R file requires no external libraries to run

#' @title Calculates and returns the LDL Value for any of the 12 different equations
#' @description This function calculates and returns the LDL value computed from any of the 12 named equations.
#' @param TC The TC (Total Cholesterol) value.
#' @param HDL The HDL (High-density lipoprotein- cholesterol) value.
#' @param TG The TG (Triglyceride) value.
#' @param EqMethod The type of equation to be used to calculate the LDL value.The type of equation to be used to calculate the LDL value.
#' EqMethod could be:("Friedewald","Ahmadi","Chen","Anandaraja","NewFormula","deCordova","Vujovic","Hattori","Puavillai","Hatta","Martin180","Martin360","Martin2000","DeLong" or "Rao").
#' @return The calculated LDL value, according to the equation of choice or a printed error message and 404,
#'  if the equation type does not exist.
#' @examples
#' LDL_eq(170.5,35.12,230,"Martin360")
#' @export
LDL_eq<-function(TC,HDL,TG,EqMethod){
  if ((EqMethod=="Friedewald")&(TG>=400)){
    message("Friedewald equation is not valid for TG values >= 400.")
    return (-1)
  }
  else if ((EqMethod=="Friedewald")&(TG<400)){
    LDL_friedwald<-TC-HDL-(TG/5)
    return(LDL_friedwald)
  }
  else if (EqMethod=="Ahmadi"){
    LDL_Ahmadi<-TC/1.19+TG/1.9-HDL/1.1-38
    return(LDL_Ahmadi)
  }
  else if (EqMethod=="Chen"){
    LDL_Chen<-(0.9*TC)-(0.9* HDL)-(0.1*TG)
    return(LDL_Chen)
  }
  else if (EqMethod=="Anandaraja"){
    LDL_Anandaraja<-(0.9*TC)-(0.9*TG/5)-28
    return(LDL_Anandaraja)
  }
  else if (EqMethod=="NewFormula"){
    LDL_NewFormula<-(0.97*TC)-(0.93*HDL)-(0.19*TG)
    return(LDL_NewFormula)
  }
  else if (EqMethod=="deCordova"){
    LDL_deCordova<-0.7516*(TC-HDL)
    return(LDL_deCordova)
  }
  else if (EqMethod=="Vujovic"){
    LDL_Vujovic<-TC-HDL-(TG/6.85)
    return(LDL_Vujovic)
  }
  else if (EqMethod=="Hattori"){
    LDL_Hattori<-(0.94*TC) - (0.94*HDL) - (0.19*TG)
    return(LDL_Hattori)
  }
  else if (EqMethod=="Puavillai"){
    LDL_Puavillai<-TC-HDL-TG/6
    return(LDL_Puavillai)
  }
  else if (EqMethod=="Hatta"){
    LDL_Hatta<- TC-HDL-TG/4
    return(LDL_Hatta)
  }
  else if (EqMethod=="Martin180"){
    LDL_Martin_180<- fLDL180(TC,HDL,TG)
    return(LDL_Martin_180)
  }
  else if (EqMethod=="Martin360"){
    LDL_Martin_360<-fLDL360(TC,HDL,TG)
    return(LDL_Martin_360)
  }
  else if (EqMethod=="Martin2000"){
    LDL_Martin_2000<-fLDL2000(TC,HDL,TG)
    return(LDL_Martin_2000)

  }
  else if (EqMethod=="DeLong"){
    LDL_DeLong<-TC-(HDL+0.16*TG)
    return(LDL_DeLong)
  }
  else if (EqMethod=="Rao"){
    LDL_Rao<-TC-HDL-(TG*(0.203-(0.00011*TG)))
    return(LDL_Rao)
  }

  else
    {
      message("No such equation type exists in the package.")
    return (404)
    }
}

#' @title Calculates and returns the LDL values using all available equations
#' @description This function calculates and returns the LDL values computed with all of the 12 named equations.
#' @param TC The TC (Total Cholesterol) value.
#' @param HDL The HDL (High-density lipoprotein-cholesterol) value.
#' @param TG The TG (Triglyceride) value.
#' @return The calculated LDL values, according to all the equations.
#' @examples
#' LDLallEq(170,35,174)
#' @export
LDLallEq <- function(TC, HDL, TG) {
  # Define the equation methods
  EqMethod <- c("Friedewald", "Ahmadi", "Chen", "Anandaraja", "NewFormula", "deCordova",
                "Vujovic", "Hattori", "Puavillai", "Hatta", "Martin180", "Martin360",
                "Martin2000", "DeLong", "Rao")

  # Initialize vectors to store results
  LDLValues <- numeric(length(EqMethod))
  MethodNames <- character(length(EqMethod))

  # Calculate LDL values for each method
  for (i in seq_along(EqMethod)) {
    LDLValues[i] <- LDL_eq(TC, HDL, TG, EqMethod[i]) # Assuming LDL_eq returns the calculated LDL value
    MethodNames[i] <- EqMethod[i]
  }

  # Create a data frame with the results
  resultsDf <- data.frame(Method = MethodNames, LDL = LDLValues)

  return(resultsDf)
}

