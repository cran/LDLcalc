# This .R file requires "utils" package to run

# Returns first value of vector
firstValue <- function(x) { return( x[length(1)] ) }

# Get TG row that corresponds to a given TG value
fGetTGRow <- function(TGvalue, TGRanges) {
  if (TGvalue<7){return (0)} # Min detected value
  if (TGvalue>13975){return (0)}# Max detected value
  TGRangesSorted = sort(c(TGvalue, TGRanges))
  index = which(TGRangesSorted==TGvalue)
  index = firstValue(index)
  index = as.integer(index/2)
  return(index)
}

# Get non-HDL column that corresponds to a given non-HDL value
fGetNonHDLcol <- function(nonHDLvalue, nonHDLRanges) {
  if (nonHDLvalue<0){return (0)}
  nonHDLRangesSorted = sort(c(nonHDLvalue, nonHDLRanges))
  index = which(nonHDLRangesSorted==nonHDLvalue)
  index = firstValue(index)
  index = as.integer(index/2)
  return(index)
}

# Get division factor for the arguments given (according to the dfCELL tables)
fGetDivFactor <- function(CHOL, HDL, TG, TGRanges, nonHDLRanges, dfCell) {
  nonHDL = CHOL - HDL
  row = fGetTGRow(TG, TGRanges)
  col = fGetNonHDLcol(nonHDL, nonHDLRanges)
  divFactor = dfCell[row, col]
  return(divFactor)
}

# Calculates the LDL value based on the Martin 180 tables and equation
fLDL180 <- function(CHOL, HDL, TG) {
  divFactor = fGetDivFactor(CHOL, HDL, TG, TGRanges180, nonHDLRanges180, dfCell180)
  LDL180 = CHOL - HDL - (TG/divFactor)
  return(round(LDL180, 0))
}

# Calculates the LDL value based on the Martin 360 tables and equation
fLDL360 <- function(CHOL, HDL, TG) {
  divFactor = fGetDivFactor(CHOL, HDL, TG, TGRanges360, nonHDLRanges360, dfCell360)
  LDL360 = CHOL - HDL - (TG/divFactor)
  return(round(LDL360, 0))
}

# Calculates the LDL value based on the Martin 2000 tables and equation
fLDL2000 <- function(CHOL, HDL, TG) {
  divFactor = fGetDivFactor(CHOL, HDL, TG, TGRanges2000, nonHDLRanges2000, dfCell2000)
  LDL2000 = CHOL - HDL - (TG/divFactor)
  return(round(LDL2000, 0))
}
