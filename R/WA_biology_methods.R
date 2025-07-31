#' @useDynLib WAFishBiology, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# WA Fish Biology methods package
# Alex Hesp Last updated August 2024
# Department of Primary Industries and Regional Development

# **************************
# General plotting functions
# **************************

#' Get maximum and minimum y values and y interval
#'
#' Used for setting plot defaults
#'
#' @keywords internal
#'
#' @param y_data y axis data for plot
#'
#' @return list object with minimum and maximum y axis values and y axis interval
Get_yaxis_scale <- function(y_data) {

  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

  ymax_data = max(1.1 * y_data)
  ymin_data = min(y_data)

  ypow = log10(ymax_data-ymin_data)
  yint = 10 ^ (ypow-round(ypow,0))

  if (yint>=0 & yint<2.5) yint = 0.2
  if (yint>=2.5 & yint<5) yint = 0.5
  if (yint>=5 & yint<7.5) yint = 1
  if (yint>=7.5) yint = 2

  yint = yint * 10^round(ypow,0) # major ticks
  ymin = yint * round(ymin_data/yint,0)
  ymax = yint * (round(ymax_data / yint,0) + 1)

  results = list(ymin = ymin,
                 ymax = ymax,
                 yint = yint)

  return(results)

}

#' Get maximum and minimum x values and x interval
#'
#' Used for setting plot defaults
#'
#' @keywords internal
#'
#' @param x_data x axis data for plot
#'
#' @return list object with minimum and maximum x axis values and x axis interval
Get_xaxis_scale <- function(x_data) {

  # modified from code provided by
  # https://peltiertech.com/calculate-nice-axis-scales-in-excel-vba/

  xmax_data = max(x_data)
  xmin_data = min(x_data)

  xpow = log10(xmax_data-xmin_data)
  xint = 10 ^ (xpow-round(xpow,0))

  if (xint>=0 & xint<2.5) xint = 0.2
  if (xint>=2.5 & xint<5) xint = 0.5
  if (xint>=5 & xint<7.5) xint = 1
  if (xint>=7.5) xint = 2

  xint = xint * 10^round(xpow,0) # major ticks
  xmin = xint * round(xmin_data/xint,0)
  xmax = xint * (round(xmax_data / xint,0) + 1)

  results = list(xmin = xmin,
                 xmax = xmax,
                 xint = xint)

  return(results)

}

#' Generic function to add axes and axes labels to plots
#'
#' @keywords internal
#' @param xmin x axis minimum
#' @param xmax x axis maximum
#' @param xint x axis tick interval
#' @param ymin y axis minimum
#' @param ymax y axis maximum
#' @param yint y axis tick interval
#' @param cexval label size
#' @param cexaxisval axis size
#' @param lwdval line widths
#' @param lineval axis offset
#' @param lasval axis orientation
#'
#' @return adds axes to plots
AddAxesAndTickLabelsToPlot <- function(xmin, xmax, xint, ymin, ymax, yint, cexval, cexaxisval, lwdval, lineval, lasval, xaxlabel, tcklen) {

  if (is.na(xmin)) xmin=0
  if (is.na(ymin)) ymin=0
  if (is.na(cexval)) cexval=1
  if (is.na(cexaxisval)) cexaxisval=1
  if (is.na(lwdval)) lwdval=1
  if (is.na(lasval)) lasval=1
  if (is.na(lineval)) lineval=0

  axis(1, at = seq(xmin, xmax, xint), line = lineval, labels = xaxlabel)
  axis(2, at = seq(ymin, ymax, yint), line = lineval, labels = F)
  axis(2, at = seq(ymin, ymax, yint), lwd=lwdval, labels=T, line=lineval, cex=cexval, cex.axis=cexaxisval, las=lasval)
}

#***********************************
# Parameter transformation functions
#***********************************


#' Inverse logit transformation
#'
#' Inverse logit transformation (back-transforms logit-transformed values, used
#' to keep parameters between 0 and 1 when fitting a model)
#'
#' @keywords internal
#'
#' @param x logit-transformed value
#'
#' @return inverse of logit-transformed value
ilogit <- function(x) {
  result = 1/(1+exp(-x))
  return(result)
}

#' Logit transformation
#'
#' Logit-transformation, used to keep parameters between 0 and 1 when fitting a model
#'
#' @keywords internal
#'
#' @param x value to be transformed
#'
#' @return logit-transformed value
logit <- function(x) {
  result = log(x/(1-x))
  return(result)
}

#**********************************
# Generic, smooth penalty functions
#**********************************

#' Logistic function (for lower bound) of penalty function for smooth upper and lower bounds (attributed to Prof. Norm Hall)
#'
#' Returns 0 if x << bound and 1 if x >> bound. The parameter, slope, determines the range over
#' which the return value changes. The function used is the logistic function, where the return
#' value is 0.5 when x == bound and 0.95 when x95 == bound + log(19)/slope
#'
#' @keywords internal
#'
#' @param x value to be check for penalty constraint
#'
#' @return   Returns 0 if x << bound and 1 if x >> bound
x_is_gt_bound <- function(x, bound, slope) {

  result = 1.0 / (1.0 + exp(slope * (bound - x)))
  return (result)
}

#' Logistic function (for upper bound) of penalty function for smooth upper and lower bounds (attributed to Prof. Norm Hall)
#'
#' Returns 0 if x << bound and 1 if x >> bound. The parameter, slope, determines the range over
#' which the return value changes. The function used is the logistic function, where the return
#' value is 0.5 when x == bound and 0.95 when x95 == bound + log(19)/slope
#'
#' @keywords internal
#'
#' @param x value to be check for penalty constraint
#'
#' @return   Returns 0 if x << bound and 1 if x >> bound
x_is_lt_bound <- function(x, bound, slope) {
  result = 1.0 - 1.0 / (1.0 + exp(slope * (bound - x)))
  return (result)
}

#' Smooth penalty function
#'
#' Smooth penalty function, equivalent functionality to ADMB posfun (excluding use of if statements)
#'
#' @keywords internal
#'
#' @param x value to be check for penalty constraint
#'
#' @return  value of x after checking for constraint
posfun <- function (x, eps, pen) {
  # Assume eps = log(19)/slope
  slope = log(19.0) / eps
  bound = eps
  result = x *  x_is_gt_bound(x, bound, slope);
  return (result)
}


#' Smooth penalty function penalty value calculation
#'
#' Smooth penalty function, equivalent functionality to ADMB posfun (excluding use of if statements)
#'
#' @keywords internal
#'
#' @param x value to be check for penalty constraint
#'
#' @return  value of x after checking for constraint
penfun <- function(x, eps, pen) {
  # The argument, pen, is not used but included for consistency with the posfun function above.
  # Assume eps = log(19)/slope
  slope = log(19.0) / eps
  bound = eps
  result = (x-eps) * (x-eps) * x_is_lt_bound(x, bound, slope)
  return (result)
}

# ***********************************************************
# Analyses associated with size-related movement growth model
# ***********************************************************

#' Calculate the probability of seeing an observed length within a juvenile habitat
#'
#' This function returns a calculated probability of seeing an observed length (x) within a juvenile habitat, for
#' use in growth model with offshore movement
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return returns probability
f1_OffMoveMod <- function(x, MeanLen, Growthsd, moveL50, moveL95) { # Opt 1 - insh

  NormProb = suppressWarnings(dnorm(x, mean = MeanLen, sd = Growthsd))
  if (is.nan(sum(NormProb))) {
    NormProb = 1e-20
  }

  ProbHab = 1 - (1 / (1 + exp(-log(19) * (x - moveL50) / (moveL95 - moveL50))))
  result = NormProb * ProbHab
  return(result)
}

#' Calculate the probability of seeing an observed length (x) within an adult  habitat
#'
#' This function returns a calculated probability of seeing an observed length within an adult habitat, for
#' use in growth model with offshore movement
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return returns probability
f2_OffMoveMod <- function(x, MeanLen, Growthsd, moveL50, moveL95) { # Opt 2 - insh

  NormProb = suppressWarnings(dnorm(x, mean = MeanLen, sd = Growthsd))
  if (is.nan(sum(NormProb))) {
    NormProb = 1e-20
  }

  ProbHab = 1 - (1 / (1 + exp(-log(19) * (x - moveL50) / (moveL95 - moveL50))))
  result = NormProb * ProbHab * x
  return(result)
}

#' Calculate the product of x.f(x), for use in calculating the expected length
#' of x, i.e. the mean, for fish in a juvenile habitat
#'
#' This function returns the product of x.f(x), for use in calculating the expected length
#' of x, i.e. the mean, for fish in a juvenile habitat
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return x.f(x)
f3_OffMoveMod <- function(x, MeanLen, Growthsd, moveL50, moveL95) { # Opt 3 - off

  NormProb = suppressWarnings(dnorm(x, mean = MeanLen, sd = Growthsd))
  if (is.nan(sum(NormProb))) {
    NormProb = 1e-20
  }

  ProbHab = 1 / (1 + exp(-log(19) * (x - moveL50) / (moveL95 - moveL50)))
  result = NormProb * ProbHab
  return(result)
}

#' Calculate the product of x.f(x), for use in calculating the expected length
#' of x, i.e. the mean, for fish in an adult habitat
#'
#' This function returns the product of x.f(x), for use in calculating the expected length
#' of x, i.e. the mean, for fish in a juvenile habitat
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return x.f(x)
f4_OffMoveMod_OffMoveMod <- function(x, MeanLen, Growthsd, moveL50, moveL95) { # Opt 4 - off

  NormProb = suppressWarnings(dnorm(x, mean = MeanLen, sd = Growthsd))
  if (is.nan(sum(NormProb))) {
    NormProb = 1e-20
  }

  ProbHab = 1 / (1 + exp(-log(19) * (x - moveL50) / (moveL95 - moveL50)))
  result = NormProb * ProbHab * x
  return(result)
}

#' Calculate the probability of a fish being in a juvenile habitat from its length
#'
#' This function calculates the probability of a fish being in a juvenile habitat from its length
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return ProbInsh
CalcProbInsh_OffMoveMod <- function(ObsLen, MeanLen, Growthsd, moveL50, moveL95) {

  a = MeanLen - (5 * Growthsd)
  b = MeanLen + (5 * Growthsd)

  # Determine the sum over all lengths of the probability
  # of seeing each length within a sample from the juvenile habitat
  temp = integrate(f1_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                   MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sum = temp$value

  # Determine the probability of seeing the observed length
  # within a sample from the juvenile habitat
  prob = f1_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)

  # Normalise the probability of seeing the observed length
  # within a sample from the juvenile habitat
  if (sum > prob) {
    ProbInsh = prob / sum
  } else {
    ProbInsh = 1E-20
  }

  return(ProbInsh)

}

#' Calculate the probability of the expected length of a fish being in an adult habitat from its length
#'
#' This function calculates the probability of a fish being in an adult habitat from its length
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return ProbOff
CalcProbOff_OffMoveMod <- function(ObsLen, MeanLen, Growthsd, moveL50, moveL95) {

  a = MeanLen - (5 * Growthsd)
  b = MeanLen + (5 * Growthsd)

  # Determine the sum over all lengths of the probability
  # of seeing each length within a sample from the adult habitat
  temp = integrate(f3_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                   MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sum = temp$value

  # Determine the probability of seeing the observed length
  # within a sample from the adult habitat
  prob = f3_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)

  # Normalise the probability of seeing the observed length
  # within a sample from the sheltered adult habitat
  ProbOff = prob / sum

  if (sum > prob) {
    ProbOff = prob / sum
  } else {
    ProbOff = 1E-20
  }

  return(ProbOff)

}

#' Calculate expected mean length of fish in the juvenile habitat, given overall mean length
#' across habitats and specified observed length
#'
#' This function calculates the expected mean length of fish in the juvenile habitat, given overall mean length
#' across habitats and specified observed length
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return MeanInsh
CalcMeanLenInsh_OffMoveMod <- function(ObsLen, MeanLen, Growthsd, moveL50, moveL95) {

  a = MeanLen - (5 * Growthsd)
  b = MeanLen + (5 * Growthsd)

  temp = integrate(f2_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                   MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sumLenProb=temp$value

  temp = integrate(f1_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                   MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sum = temp$value
  MeanInsh = sumLenProb / sum

  return(MeanInsh)

}

#' Calculate expected mean length of fish in the adult habitat, given overall mean length
#' across habitats and specified observed length
#'
#' This function calculates the expected mean length of fish in the adult habitat, given overall mean length
#' across habitats and specified observed length
#'
#' @keywords internal
#'
#' @param params x, MeanLen, Growthsd, moveL50, moveL95
#'
#' @return MeanOff
CalcMeanLenOff_OffMoveMod <- function(ObsLen, MeanLen, Growthsd, moveL50, moveL95) {

  a = MeanLen - (5 * Growthsd)
  b = MeanLen + (5 * Growthsd)

  temp = integrate(f4_OffMoveMod_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                   MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sumLenProb=temp$value

  temp <- integrate(f3_OffMoveMod, lower = a, upper = b, rel.tol = 1E-12,
                    MeanLen = MeanLen, Growthsd = Growthsd, moveL50 = moveL50, moveL95 = moveL95)
  sum = temp$value

  MeanOff = sumLenProb / sum

  return(MeanOff)

}

#' Calculate negative log-likelihood associated with length-at-age data, given
#' parameters for size-related movement growth model
#'
#' This function negative log-likelihood associated with length-at-age data, given
#' parameters for size-related movement growth model (as described by Hesp et al., 2004)
#'
#' @keywords internal
#'
#' @param params
#'
#' @return NLL
CalcNLL_OffMoveMod <- function(params) {

  # get params
  lnmoveL50 = params[1]
  lnmoveSlope = params[2]
  lnLinf = params[3]
  lnvbK = params[4]
  tzero = params[5]
  CV = exp(params[6])

  # parameter penalties
  eps = 0.001
  param_pen=0
  L50MovePen=0
  SlopeMovePen=0
  Linf_pen=0
  vbK_pen=0
  tzero_pen=0
  CV_pen = 0

  # L50Move
  L50MovePen = penfun(lnmoveL50 - lnL50MoveLw, eps, L50MovePen);
  temp_parm = lnL50MoveLw + posfun(lnmoveL50 - lnL50MoveLw, eps, L50MovePen);
  L50MovePen = L50MovePen + penfun(lnL50MoveUp - temp_parm, eps, L50MovePen);
  temp_parm = lnL50MoveUp - posfun(lnL50MoveUp - temp_parm, eps, L50MovePen);
  moveL50 = exp(temp_parm);
  param_pen = param_pen + 1000*L50MovePen;

  # SlopeMove
  SlopeMovePen = penfun(lnmoveSlope - lnSlopeMoveLw, eps, SlopeMovePen);
  temp_parm = lnSlopeMoveLw + posfun(lnmoveSlope - lnSlopeMoveLw, eps, SlopeMovePen);
  SlopeMovePen = SlopeMovePen + penfun(lnSlopeMoveUp - temp_parm, eps, SlopeMovePen);
  temp_parm = lnSlopeMoveUp - posfun(lnSlopeMoveUp - temp_parm, eps, SlopeMovePen);
  moveSlope = exp(temp_parm);
  param_pen = param_pen + 1000*SlopeMovePen
  moveL95 <<- (log(19) / moveSlope) + moveL50

  # Linf
  Linf_pen = penfun(lnLinf - lnLinfLw, eps, Linf_pen);
  temp_parm = lnLinfLw + posfun(lnLinf - lnLinfLw, eps, Linf_pen);
  Linf_pen = Linf_pen + penfun(lnLinfUp - temp_parm, eps, Linf_pen);
  temp_parm = lnLinfUp - posfun(lnLinfUp - temp_parm, eps, Linf_pen);
  Linf = exp(temp_parm);
  param_pen = param_pen + 1000*Linf_pen;

  # vbK
  vbK_pen = penfun(lnvbK - lnvbKLw, eps, vbK_pen);
  temp_parm = lnvbKLw + posfun(lnvbK - lnvbKLw, eps, vbK_pen);
  vbK_pen = vbK_pen + penfun(lnvbKUp - temp_parm, eps, vbK_pen);
  temp_parm = lnvbKUp - posfun(lnvbKUp - temp_parm, eps, vbK_pen);
  vbK = exp(temp_parm);
  param_pen = param_pen + 1000*vbK_pen;

  # tzero
  tzero_pen = penfun(tzero - tzeroLw, eps, tzero_pen);
  temp_parm = tzeroLw + posfun(tzero - tzeroLw, eps, tzero_pen);
  tzero_pen = tzero_pen + penfun(tzeroUp - temp_parm, eps, tzero_pen);
  temp_parm = tzeroUp - posfun(tzeroUp - temp_parm, eps, tzero_pen);
  tzero = temp_parm;
  param_pen = param_pen + 1000*tzero_pen;

  # cv
  CV_pen = penfun(CV - CVLw, eps, CV_pen);
  temp_parm = CVLw + posfun(CV - CVLw, eps, CV_pen);
  CV_pen = CV_pen + penfun(CVUp - temp_parm, eps, CV_pen);
  temp_parm = CVUp - posfun(CVUp - temp_parm, eps, CV_pen);
  CV = temp_parm;
  param_pen = param_pen + 1000*CV_pen;

  ExpLen = rep(0,SampleSize)
  NLL=0
  for (i in 1:SampleSize) {

    ObsAge = ObsAges[i]
    ObsLen = ObsLengths[i]

    # est mean length at age from vb growth equation
    MeanLen = Linf * (1-exp(-vbK * (ObsAge - tzero)))

    # get sd, for mean length at age
    Growthsd = CV * MeanLen

    # est mean length of each fish at age, given habitat
    if (Habitat[i] == 0) {
      ExpLen[i] = CalcMeanLenInsh_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
      ObsProbInsh = CalcProbInsh_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
      NLL1 = -log(ObsProbInsh)
    } else {
      ExpLen[i] = CalcMeanLenOff_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
      ObsProbOff = CalcProbOff_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
      NLL1 = -log(ObsProbOff)
    }

    NLL = NLL + NLL1 + param_pen

  }
  cat("NLL",NLL,"Linf",Linf,"vbK",vbK,"tzero",tzero,"CV",CV,"L50",moveL50,"Slope",moveSlope, "Pen",param_pen,'\n')

  return(NLL)

}

#' Fits growth model allowing for size-related movements of fish between 2 habitats
#'
#' This function growth model allowing for size-related movements of fish between 2 habitats,
#' (see Hesp et al., 2004)
#'
#' @param params input param values
#' @param SampleSize observed sample size for length/age/habitat data
#' @param ObsAges observed ages
#' @param ObsLengths observed lengths
#' @param Habitat observed habitat 0-first habitat, 1=final habitat
#' @param lnL50MoveUp upper bound movement parameter
#' @param lnL50MoveLw lower bound movement parameter
#' @param lnSlopeMoveUp upper bound movement parameter
#' @param lnSlopeMoveLw lower bound movement parameter
#' @param lnLinfUp upper bound growth parameter
#' @param lnLinfLw lower bound growth parameter
#' @param lnvbKUp upper bound growth parameter
#' @param lnvbKLw lower bound growth parameter
#' @param tzeroUp upper bound growth parameter
#' @param tzeroLw lower bound growth parameter
#' @param CVLw upper bound growth variation parameter
#' @param CVUp lower bound growth variation parameter
#'
#' @return nll, convergence, SampleSize, ParamEst, params, vcov.params, cor.params, IndivFishRes
#' @examples
#' library(L3Assess)
#' library(WAFishBiology)
#' # simulate some standard length-at-age data
#' set.seed(123)
#' SampleSize=1000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 20
#' TimeStep = 1/12 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/12
#' FishMort = 0.1
#' MaxLen = 800
#' LenInc = 10
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityAtLen = NA # selectivity vector
#' SelParams = c(100, 20) # L50, L95-L50 for gear selectivity
#' RetenParams = c(NA, NA) # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 300
#' vbK = 0.5
#' CVSizeAtAge = 0.03
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                                SelParams, RetenParams, SelectivityAtLen, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#'
#' ObsAges = Res$ObsDecAgeRetCatch
#' ObsLengths = Res$ObsRandLenRetCatch
#' plot(ObsAges, ObsLengths, xlim=c(0,20), ylim=c(0,400))
#' # Assign each fish to inshore or offshore habitat
#' moveL50 = 200
#' moveSlope = 0.05
#' moveL95 = (log(19)/moveSlope) + moveL50
#' Habitat <- rep(0,SampleSize)
#' for (i in 1:SampleSize) {
#'   randnum = runif(1,0,1)
#'   tempProbOff = 1 / (1 + exp(-log(19) * (ObsLengths[i] - moveL50) / (moveL95 - moveL50)))
#'   if (randnum<tempProbOff) {
#'     Habitat[i] = 1
#'   } else {
#'     Habitat[i] = 0
#'   }
#' }
#' # Reduce relative sampling intensity in one area (selectivity offshore)
#' # at is unlikely that sampling will be of the same intensity in both habitats
#' x = which(Habitat==0)
#' y = which(Habitat==1)
#' RandAgeInsh = ObsAges[x]
#' RandAgeOff = ObsAges[y]
#' RandLenInsh = ObsLengths[x]
#' RandLenOff = ObsLengths[y]
#' y=which(Habitat==1)
#' length(y)
#' randnum=runif(length(y),0,1)
#' yy=which(randnum<0.25)
#' RandAgeOff = RandAgeOff[yy]
#' RandLenOff = RandLenOff[yy]
#' par(mfrow=c(2,2),mar=c(4,3,2,2))
#' plot(RandAgeInsh,RandLenInsh, xlim=c(0,20), ylim=c(0,500), col="red")
#' plot(RandAgeOff,RandLenOff, xlim=c(0,20), ylim=c(0,500), col="blue")
#' plot(RandAgeOff,RandLenOff, xlim=c(0,20), ylim=c(0,500), col="blue")
#' points(RandAgeInsh,RandLenInsh, col="red")
#' plot(RandAgeInsh,RandLenInsh, xlim=c(0,20), ylim=c(0,500), col="red")
#' points(RandAgeOff,RandLenOff, col="blue")
#' # Insh
#' lbnd=Res$lbnd
#' InshLenCat = trunc(RandLenInsh/10)*10
#' histdat = hist(InshLenCat, breaks=c(lbnd,800), plot=F)
#' InshLenFreq = histdat$counts
#' breaks = histdat$breaks
#' plot(lbnd,InshLenFreq,"l")
#' # Off
#' OffLenCat = trunc(RandLenOff/10)*10
#' histdat = hist(OffLenCat, breaks=c(lbnd,800), plot=F)
#' OffLenFreq = histdat$counts
#' breaks = histdat$breaks
#' lines(lbnd,OffLenFreq,"l",col="blue")
#' # plot lengths at specified ages. Lengths of fish in offshore habitat
#' # expected to be larger at younger ages
#' SpecAge = 1
#' lbnd=Res$lbnd
#' x=which(ObsAges >= SpecAge & ObsAges < SpecAge+1 & Habitat == 0)
#' InshLenCat = trunc(ObsLengths[x]/10)*10
#' histdat = hist(InshLenCat, breaks=c(lbnd,800), plot=F)
#' InshLenFreq = histdat$counts
#' breaks = histdat$breaks
#' x=which(ObsAges >= SpecAge & ObsAges < SpecAge+1 & Habitat == 1)
#' OffLenCat = trunc(ObsLengths[x]/10)*10
#' histdat = hist(OffLenCat, breaks=c(lbnd,800), plot=F)
#' OffLenFreq = histdat$counts
#' breaks = histdat$breaks
#' ymax=1.2*max(c(InshLenFreq, OffLenFreq))
#' plot(lbnd,InshLenFreq,"l",ylim=c(0,ymax))
#' lines(lbnd,OffLenFreq, col="blue")
#' # Fit size-related movement growth model
#' # set parameter bounds
#' lnL50MoveUp = log(300)
#' lnL50MoveLw = log(50)
#' lnSlopeMoveUp = log(0.5)
#' lnSlopeMoveLw = log(0.01)
#' lnLinfUp = log(350)
#' lnLinfLw = log(250)
#' lnvbKUp = log(0.8)
#' lnvbKLw = log(0.2)
#' tzeroLw = -5
#' tzeroUp = 1
#' CVLw = 0.02
#' CVUp = 0.15
#' # set starting values for parameters
#' InitmoveL50 = rnorm(1,200,20); InitmoveL50
#' InitmoveSlope = rnorm(1,0.05,0.002); moveSlope
#' InitLinf = rnorm(1,300,20); InitLinf
#' InitvbK = rnorm(1,0.5,0.05); InitvbK
#' Inittzero = rnorm(1,0,0.05); Inittzero
#' InitGrowthcv = rnorm(1,0.1,0.005); InitGrowthcv
#' params = c(log(InitmoveL50), log(InitmoveSlope),log(InitLinf), log(InitvbK), Inittzero,log(InitGrowthcv))
#' # fit model (takes a while, uses nlminb and Amoeba routine for optimisation)
#' FittedRes=GetOffMoveGrowthModResults(params, SampleSize, ObsAges, ObsLengths, Habitat, lnL50MoveUp, lnL50MoveLw,
#'                                      lnSlopeMoveUp, lnSlopeMoveLw, lnLinfUp, lnLinfLw, lnvbKUp, lnvbKLw,
#'                                      tzeroLw, tzeroUp, CVLw, CVUp)
#' # plot expected mean lengths at age for each habitat
#' par(mfrow=c(2,2),mar=c(5,4,2,2))
#' x=which(FittedRes$IndivFishRes$Habitat==0)
#' plot(FittedRes$IndivFishRes$ObsAges[x],FittedRes$IndivFishRes$ObsLengths[x], xlim=c(0,20), ylim=c(0,400),
#'      col="red",xlab="Age",ylab="Length",bty='n', cex.main=0.8, main="Juv Habitat")
#' points(FittedRes$IndivFishRes$ObsAges[x], FittedRes$IndivFishRes$ExpLen[x], col="black")
#'
#' y=which(FittedRes$IndivFishRes$Habitat==1)
#' plot(FittedRes$IndivFishRes$ObsAges[y],FittedRes$IndivFishRes$ObsLengths[y], xlim=c(0,20), ylim=c(0,400),
#'      col="blue",xlab="Age",ylab="Length",bty='n', cex.main=0.8, main="Adult Habitat")
#' points(FittedRes$IndivFishRes$ObsAges[y], FittedRes$IndivFishRes$ExpLen[y], col="black")
#'
#' plot(FittedRes$IndivFishRes$ObsAges[x], FittedRes$IndivFishRes$ExpLen[x], xlim=c(0,20), ylim=c(0,400),
#'      xlab="Age",ylab="Length",bty='n', col="red")
#' points(FittedRes$IndivFishRes$ObsAges[y], FittedRes$IndivFishRes$ExpLen[y], col="blue")
#' @export
GetOffMoveGrowthModResults <- function(params, SampleSize, ObsAges, ObsLengths, Habitat, lnL50MoveUp, lnL50MoveLw,
                                       lnSlopeMoveUp, lnSlopeMoveLw, lnLinfUp, lnLinfLw, lnvbKUp, lnvbKLw,
                                       tzeroLw, tzeroUp, CVLw, CVUp) {

  # fit model with nlminb
  cat("fitting with nlminb",'\n')
  nlmb <- nlminb(params, CalcNLL_OffMoveMod, gradient = NULL, hessian = TRUE)
  nlmb$objective

  # Run Nelder-Mead optimization
  params <- nlmb$par
  cat("fitting with Optim: Nelder-Mead",'\n')
  Amoeb <- optim(params, CalcNLL_OffMoveMod, method = "Nelder-Mead")
  params <- Amoeb$par

  # get estimates
  Amoeb$value # value of nll
  Amoeb$convergence
  Amoeb$par

  # calculate uncertainty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(Amoeb$par, CalcNLL_OffMoveMod)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  EstmoveL50 = c(exp(Amoeb$par[1]), exp(Amoeb$par[1] + c(-1.96, 1.96) * ses[1]))
  EstmoveSlope = c(exp(Amoeb$par[2]), exp(Amoeb$par[2] + c(-1.96, 1.96) * ses[2]))
  EstmoveL95 = (log(19) / EstmoveSlope) + EstmoveL50
  EstLinf <- c(exp(Amoeb$par[3]), exp(Amoeb$par[3] + c(-1.96, 1.96) * ses[3]))
  EstvbK <- c(exp(Amoeb$par[4]), exp(Amoeb$par[4] + c(-1.96, 1.96) * ses[4]))
  Esttzero <- c(Amoeb$par[5], Amoeb$par[5] + c(-1.96, 1.96) * ses[5])
  EstCV <- c(exp(Amoeb$par[6]), exp(Amoeb$par[6] + c(-1.96, 1.96) * ses[6]))

  ParamEst = t(data.frame(moveL50=round(EstmoveL50,1), moveSlope=round(EstmoveSlope,3),
                          DerivedL95=round(EstmoveL95,1),
                          Linf=round(EstLinf,1), vbK=round(EstvbK,3),
                          tzero=round(Esttzero,2), CV=round(EstCV,3)))
  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")


  # store value of objective function
  nll = Amoeb$value

  # store convergence value
  convergence = Amoeb$convergence

  # get expected length for each fish, in its habitat
  ExpLen = rep(NA,SampleSize)
  for (i in 1:SampleSize) {

    ObsAge = ObsAges[i]
    ObsLen = ObsLengths[i]

    # est mean length at age from vb growth equation
    Linf = ParamEst[4,1]
    vbK = ParamEst[5,1]
    tzero = ParamEst[6,1]

    MeanLen = Linf * (1-exp(-vbK * (ObsAge - tzero)))

    # get sd, for mean length at age
    CV = ParamEst[7,1]
    Growthsd = CV * MeanLen
    moveL50 = ParamEst[1,1]
    moveL95 = ParamEst[3,1]

    # est mean length of each fish at age, given habitat
    if (Habitat[i] == 0) {
      ExpLen[i] = CalcMeanLenInsh_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
    } else {
      ExpLen[i] = CalcMeanLenOff_OffMoveMod(ObsLen, MeanLen, Growthsd, moveL50, moveL95)
    }
  }

  IndivFishRes = data.frame(ObsAges=ObsAges, ObsLengths=ObsLengths,
                            Habitat=Habitat, ExpLen=ExpLen)

  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.params = vcov.params,
                 cor.params = cor.params,
                 IndivFishRes = IndivFishRes)


  return(results)

}


# *********************************************************
# Analyses of length effect on spawning duration and timing
# *********************************************************

#' Calculate daily proportions of spawning capable fish of a given length, using spawning duration model
#'
#' This function returns calculated values for the daily proportions of spawning capable fish
#' of a specified length, using spawning duration model (i.e. asymmetric double logistic model with
#' variable height)
#'
#' @param FishLen specified fish length
#' @param DistnType DecDay specified decimal days of the year
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#'
#' @return returns a list object containing expected proportions of spawning capable fish, of a given length,
#' for each specified day of the year, inputted as decimal day (P_t_s), and the value of a penalty
#' function used to ensure the maximum proportion of spawning capable fishing is always between 0 and 1
#' @examples
#' FishLen = 500
#' DecDay = seq(0,1,0.01)
#' lnL50 = log(300)
#' lnslope= log(0.015)
#' logitPkSpawn = log(0.6/(1-0.6))
#' lnkappa = log(0.0002)
#' lnkappa2 = log(0.0002)
#' lnslope1 = log(30)
#' lnslope2 = log(30)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' res=CalcDailySpawnProps_SpDurMod(FishLen, DecDay, params)
#' plot(DecDay,res$P_t_s,"l")
#' @export
CalcDailySpawnProps_SpDurMod <- function(fish_size, DecDay, params) {

  param_pen = 0
  temp_parm = exp(params[2])
  eps = 0.001
  param_pen = 0
  uppbound = 0.5 # equals slope of about 0.5
  param_pen = param_pen + penfun(uppbound - temp_parm, eps)
  temp_parm = uppbound - posfun(uppbound - temp_parm, eps, param_pen)
  # cat("temp_parm",temp_parm,"uppbound",uppbound,"param_pen",param_pen,'\n')

  L50 = exp(params[1])
  slope = exp(params[2])
  Peak_spawn =  1/(1+exp(-params[3])) # ilogit transform
  kappa = exp(params[4])
  kappa2 = exp(params[5])
  slope1 = exp(params[6])
  slope2 = exp(params[7])

  d50 = Peak_spawn - kappa * fish_size
  d50_2 = Peak_spawn + kappa2 * fish_size
  Height = 1 / (1 + exp(-slope*(fish_size-L50)))

  P_t_s = 1 / (1 + exp(-slope1*(DecDay-d50))) *
    1 / (1 + exp(slope2*(DecDay-d50_2))) * Height

  results = list(P_t_s=P_t_s,
             param_pen=param_pen,
             d50=d50,
             d50_2=d50_2,
             Height=Height)

  return(results)

}

#' Calculate daily proportions of spawning capable fish of a given length, using spawning duration model,
#' used for calculating NLL
#'
#' This function returns calculated values for the daily proportions of spawning capable fish
#' of a specified length, using spawning duration model (i.e. asymmetric double logistic model with
#' variable height)
#'
#' @keywords internal
#'
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#'
#' @return returns a list object containing expected proportions of spawning capable fish, of a given length,
#' for each specified day of the year, inputted as decimal day (P_t_s), and the value of a penalty
#' function used to ensure the maximum proportion of spawning capable fishing is always between 0 and 1
CalcDailySpawnProps2_SpDurMod <- function(params) {

  param_pen = 0
  temp_parm = exp(params[2])
  eps = 0.001
  param_pen = 0
  uppbound = 0.5 # equals slope of about 0.5
  param_pen = param_pen + penfun(uppbound - temp_parm, eps)
  temp_parm = uppbound - posfun(uppbound - temp_parm, eps, param_pen)
  # cat("temp_parm",temp_parm,"uppbound",uppbound,"param_pen",param_pen,'\n')

  L50 = exp(params[1])
  slope = exp(params[2])
  Peak_spawn =  1/(1+exp(-params[3])) # ilogit transform
  kappa = exp(params[4])
  kappa2 = exp(params[5])
  slope1 = exp(params[6])
  slope2 = exp(params[7])

  d50 = Peak_spawn - kappa * FishLen
  d50_2 = Peak_spawn + kappa2 * FishLen
  Height = 1 / (1 + exp(-slope*(FishLen-L50)))

  P_t_s = 1 / (1 + exp(-slope1*(DecDay-d50))) *
    1 / (1 + exp(slope2*(DecDay-d50_2))) * Height

  results = list(P_t_s=P_t_s,
             param_pen=param_pen)

  return(results)

}

#' Simulated data with daily proportions of spawning capable fish for a species, based on using spawning duration model
#'
#' This function simulates data containing daily proportions of spawning capable fish for a specified with
#' specified growth characteristics, and spawning characteristics according to values of specified parameters
#' for the spawning duration model (i.e. asymmetric double logistic model with variable height)
#'
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#' @param GrowthEqn 1=von Bertalanffy, 2=Schnute, 3=Somers seasonal curve
#' @param nSamples number of samples
#' @param nSexes number of sexes 1=single or combined sex, 2=separate sexes
#' @param MinAge  minimum age
#' @param MaxAge maximum age
#' @param AgeStep age interval
#' @param Linf asymptotic length
#' @param vbK von Bertalanffy growth coefficient
#' @param tzero hypothetical length at age zero
#' @param Growth_sd standard deviation of lengths at age
#'
#' @return returns a list object containing input data for dates of fish capture in decimal days (DecDay),
#' month of capture (MM), corresponding fish lengths (FishLen), expected probability of each fish spawning
#' according to its date of capture and size (ExpSpawnProbs), and randomly-generated maturity status (ObsMatStatus, 0=immature,
#' 1=mature)
#' @examples
#' # Simulate spawning proportion data
#' # First, simulate length and age data
#' set.seed(123)
#' GrowthEqn=1 # von Bertalanffy
#' nSamples = 5000
#' nSexes = 1
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = NA
#' Linf = 1000
#' vbK = 0.1
#' tzero = 0
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = 0.08
#' # Spawning duration model parameters
#' lnL50 = log(300)
#' lnslope= log(0.015)
#' logitPkSpawn = log(0.6/(1-0.6))
#' lnkappa = log(0.0002)
#' lnkappa2 = log(0.0002)
#' lnslope1 = log(30)
#' lnslope2 = log(30)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' Res=SimulateSpawningDurationData(params, GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
#'                                  AgeStep, Ref_ages, Growth_params, Growth_cv)
#' @export
SimulateSpawningDurationData <- function(params, GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
                                         AgeStep, Ref_ages, Growth_params, Growth_cv) {

  Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
                                AgeStep, Ref_ages, Growth_params, Growth_cv)

  # Get expected probability of spawning, given size of fish and time of year
  FishLen=round(Res$ObsLen,0)
  DecDay = Res$ObsAge - trunc(Res$ObsAge/1)*1
  MM = 1+((trunc(DecDay/(1/12))*(1/12))*12)

  res=CalcDailySpawnProps_SpDurMod(FishLen, DecDay, params)
  ExpSpawnProbs=res$P_t_s

  # assign mature /  non-mature, based on probability and random numbers
  ObsMatStatus = rep(NA,nSamples)
  for (i in 1:nSamples) {
    randnum = runif(1,0,1) # draw random number, between 0 and 1, from uniform distribution
    if (randnum <= ExpSpawnProbs[i]) {
      ObsMatStatus[i] = 1
    } else {
      ObsMatStatus[i] = 0
    }
  }

  results = list(FishLen=FishLen,
                 DecDay=DecDay,
                 MM=MM,
                 ExpSpawnProbs=ExpSpawnProbs,
                 ObsMatStatus=ObsMatStatus)

  return(results)
}

#' Calculate observed monthly proportions of spawning capable fish
#'
#' This function calculates monthly proportions of spawning capable fish for observed
#' input data on months and observed maturity status
#'
#' @keywords internal
#'
#' @param subDat data frame with observed months (MM) and maturity status (ObsMatStatus,
#' 0=immature, 1=mature)
#'
#' @return observed monthly proportions of spawning capable fish (Probs)
CalcMonthlyObsSpawnProps_SpDurMod <- function(subDat) {

  Probs = rep(0, length(seq(1,12,1))) # calculate monthly
  for (j in seq(1,12,1)) {
    xx = which(subDat$MM >= j-0.0001 & subDat$MM <= j+0.0001)
    Probs[j] = length(which(subDat$ObsMatStatus[xx] == 1)) / length(xx)
    # cat("j",j,"n",length(xx),"probs",Probs[j],'\n')
  }
  return(Probs)
}

#' Plot daily proportions of spawning capable fish for specified length categories
#'
#' This function plots calculated values for the daily proportions of spawning capable fish
#' for specified length categories and estimated daily proportions from spawning duration model,
#' (i.e. asymmetric double logistic model with variable height), if parameters are specified
#'
#' @param ObsSpawnDat list including observed fish lengths (FishLen), observed months of capture (MM)
#' and observed maturity status categories (ObsMatStatus, 0=immature, 1=mature)
#' @param lbnds lower bounds for length classes to be plotted
#' @param ubnds upper bounds for length classes to be plotted
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#'
#' @return plots
#' @examples
#' # Simulate spawning proportion data
#' # First, simulate length and age data
#' set.seed(123)
#' GrowthEqn=1 # von Bertalanffy
#' nSamples = 5000
#' nSexes = 1
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = NA
#' Linf = 1000
#' vbK = 0.1
#' tzero = 0
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = 0.08
#' # Spawning duration model parameters
#' lnL50 = log(300)
#' lnslope= log(0.015)
#' logitPkSpawn = log(0.6/(1-0.6))
#' lnkappa = log(0.0002)
#' lnkappa2 = log(0.0002)
#' lnslope1 = log(30)
#' lnslope2 = log(30)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' Res=SimulateSpawningDurationData(params, GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
#'                                  AgeStep, Ref_ages, Growth_params, Growth_cv)
#' ObsMatStatus=Res$ObsMatStatus
#' FishLen=Res$FishLen
#' MM=Res$MM
#' ObsSpawnDat = data.frame(MM=MM,FishLen=FishLen,ObsMatStatus=ObsMatStatus)
#' # params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' params = NA
#' lbnds = seq(300,900,100)
#' ubnds = lbnds + 100
#' par(mfrow=c(3,2), mar=c(4,2,2,2))
#' PlotSpawningDurationData(ObsSpawnDat, lbnds, ubnds, params)
#' @export
#'
PlotSpawningDurationData <- function(ObsSpawnDat, lbnds, ubnds, params) {

  nLenCats = length(lbnds)
  MMabb = substr(month.abb, 1, 1)
  for (i in 1:nLenCats) {
    subDat = ObsSpawnDat[ObsSpawnDat$FishLen >= lbnds[i] & ObsSpawnDat$FishLen < ubnds[i],]
    Probs = CalcMonthlyObsSpawnProps_SpDurMod(subDat)
    plot(seq(0.5,11.5,1), Probs, xaxt='n', yaxt="n", cex=0.8, cex.axis=1,bty='n', ylim=c(0,1), xlim=c(0,12),
         xlab = list(" Month",cex=1.2), ylab=list("Prop. spawning",cex=1.2), main=paste0(lbnds[i],"-",ubnds[i]," mm"), cex.main=1)
    AddAxesAndTickLabelsToPlot(xmin=0.5, xmax=11.5, xint=1, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0, lasval=1, xaxlabel = MMabb[1:12], tcklen = 0.03)
    if(!is.na(params[i])) {
      FishLen=(lbnds[i] + ubnds[i])/2
      DecDay_plot = seq(0,1,0.01)
      Res=CalcDailySpawnProps_SpDurMod(FishLen, DecDay_plot, params)
      lines(DecDay_plot*12,Res$P_t_s,col=i)
    }
  } # i
}

#' Calculate NLL for spawning duration model given data and parameter values
#'
#' Calculate negative log-likelihood for spawning duration model given data and parameter values
#'
#' @keywords internal
#'
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#'
#' @return negative log-likelihood (NLL)
CalcNLL_SpDurMod <- function(params) {

  Res = CalcDailySpawnProps2_SpDurMod(params)
  Likelihood = rep(-999, length(DecDay))
  Likelihood[which(ObsMatStatus==1)] = Res$P_t_s[which(ObsMatStatus==1)]
  Likelihood[which(ObsMatStatus==0)] = 1 - Res$P_t_s[which(ObsMatStatus==0)]
  LL <- log(Likelihood + 1E-4)
  NLL = -sum(LL) + Res$param_pen
  results = NLL
  cat("NLL",NLL,"params",params,"Res$param_pen",Res$param_pen,'\n')

  return(results)

}

#' Get results for fitted spawning duration model
#'
#' This function fits a model for estimating spawning season timing and duration, allowing for
#' fish length effects. The curves associated with the fitted model, for fish of a given length,
#' have the appearance are approximately bell-shaped, and can vary in heights among fish of different
#' lengths. The model may also described as an asymmetric double logistic model with variable height.
#'
#' @param params specified (or estimated), transformed parameter values, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#' @param FishLen Fish lengths
#' @param DecDay Days of year of fish capture (in decimal years)
#' @param ObsMatStatus maturity status, 0=immature, 1=mature
#' @param MinLen Minimum fish length (for determining lengths to output expected model curves)
#' @param MaxLen Maximum fish length for analysis
#' @param LenInc Length increment (for determining lengths to output expected model curves)
#' @param nsims number of samples to generate using parametric resampling. Set to zero to turn off resampling.
#'
#' @return estimated parameters (params), negative log-likelihood of fitted model (NLL),
#' convergence statistic (0=converged), estimated variance-covariance matrix for estimated parameters,
#' estimated parameters with 95 percent asymptotic confidence limits (including back-transformed parameters),
#' estimated parameters with 95 percent confidence limits dervied from parametric resampling of estimated
#' parameters and associated variance-covariance matrix (including back-transformed parameters),
#'
#' @examples
#' # Simulate spawning proportion data
#' set.seed(123)
#' GrowthEqn=1 # von Bertalanffy
#' nSamples = 5000
#' nSexes = 1
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = NA
#' Linf = 1000
#' vbK = 0.1
#' tzero = 0
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = 0.08
#' # Spawning duration model parameters
#' lnL50 = log(300)
#' lnslope= log(0.015)
#' logitPkSpawn = log(0.6/(1-0.6))
#' lnkappa = log(0.0002)
#' lnkappa2 = log(0.0002)
#' lnslope1 = log(30)
#' lnslope2 = log(30)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' Res=SimulateSpawningDurationData(params, GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
#'                                  AgeStep, Ref_ages, Growth_params, Growth_cv)
#' ObsMatStatus=Res$ObsMatStatus
#' FishLen=Res$FishLen
#' DecDay=Res$DecDay
#' # Spawning duration model parameters - initial values
#' lnL50 = log(350)
#' lnslope= log(0.02)
#' logitPkSpawn = log(0.7/(1-0.7))
#' lnkappa = log(0.00015)
#' lnkappa2 = log(0.00035)
#' lnslope1 = log(25)
#' lnslope2 = log(35)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' nsims=20
#' MinLen=0
#' MaxLen=1100
#' LenInc=100
#' res=GetSpawningDurationModelResults(params, FishLen, DecDay, ObsMatStatus, MinLen, MaxLen, LenInc, nsims)
#' @export
GetSpawningDurationModelResults <- function(params, FishLen, DecDay, ObsMatStatus, MinLen, MaxLen, LenInc, nsims) {

  FishLen=FishLen
  DecDay=DecDay
  ObsMatStatus=ObsMatStatus

  nlmb <- nlminb(params, CalcNLL_SpDurMod, gradient = NULL, hessian = TRUE,
                 control=list(trace=1, rel.tol=0.0000001))

  params = nlmb$par
  NLL = nlmb$objective
  convergence = nlmb$convergence

  # get variance-covariance matrix, from fitted model
  hess.out = optimHess(params, CalcNLL_SpDurMod)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params))
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  # get parameter estimates and asymptotic error estimates
  EstL50 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
  Estslope = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
  temp=c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
  EstPkSpawn = exp(temp) / (exp(temp)+1) # ilogit transform
  Estkappa = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
  Estkappa2 = c(exp(nlmb$par[5]), exp(nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
  Estslope1 = c(exp(nlmb$par[6]), exp(nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
  Estslope2 = c(exp(nlmb$par[7]), exp(nlmb$par[7] + c(-1.96, 1.96) * ses[7]))

  ParamEst = t(data.frame(EstL50=round(EstL50,2), Estslope=round(Estslope,3),
                          EstPkSpawn=round(EstPkSpawn,3), Estkappa=round(Estkappa,5),
                          Estkappa2=round(Estkappa2,5), Estslope1=round(Estslope1,3),
                          Estslope2=round(Estslope2,3)))
  ParamEst.sim = NA
  if (nsims > 0) {
    # Use resampling approach to get estimates of uncertainty of derived outputs for model
    sims = data.frame(MASS::mvrnorm(n = nsims, params, vcov.params))
    names(sims) = c("ln_L50", "ln_slope","lgt_PkSpawn","ln_kappa","ln_kappa2",
                    "ln_slope1","ln_slope2")

    sims$L50 = exp(sims$ln_L50)
    sims$slope = exp(sims$ln_slope)
    sims$Peak_spawn = exp(sims$lgt_PkSpawn) / (exp(sims$lgt_PkSpawn)+1) # ilogit transform
    sims$kappa = exp(sims$ln_kappa)
    sims$kappa2 = exp(sims$ln_kappa2)
    sims$slope1 = exp(sims$ln_slope1)
    sims$slope2 = exp(sims$ln_slope2)

    # Recalculate 95% CLs for estimated parameters and back-transformed values of estimated parameters
    sims.mean = apply(sims[, 1:14], MARGIN=2, function(x) mean(x))
    sims.median = apply(sims[, 1:14], MARGIN=2, function(x) quantile(x, 0.5))
    sims.lowCL = apply(sims[, 1:14], MARGIN=2, function(x) quantile(x, 0.025))
    sims.uppCL = apply(sims[, 1:14], MARGIN=2, function(x) quantile(x, 0.975))
    ParamEst.sim = round(cbind(sims.mean, sims.median, sims.lowCL, sims.uppCL), 3)
  }

  # get relationships between maximum proportion spawning and spawning duration with fish length,
  # with uncertainty
  SizeIntPlot = MaxLen / 100 # calculate for 100 sizes, from 0 to maxlen
  FishLen = seq(0,MaxLen,SizeIntPlot)
  nFishLen = length(FishLen)
  nDecDay = length(DecDay)
  P_t_s.sim <- array(dim=c(nsims,nFishLen,nDecDay))
  MaxCurveHeightEst = data.frame(matrix(nrow=nsims,ncol=nFishLen))
  colnames(MaxCurveHeightEst) = FishLen
  MaxCurveHeightEst = as.matrix(MaxCurveHeightEst)
  SpawnDurEst = MaxCurveHeightEst
  MaxPropSpawnVsFishLen=NA
  SpawnDurVsFishLen=NA
  if (nsims > 0) {
    for (i in 1:nsims) {

      d50 = sims$Peak_spawn[i] - sims$kappa[i] * FishLen
      d50_2 = sims$Peak_spawn[i] + sims$kappa2[i] * FishLen
      Height = 1 / (1 + exp(-sims$slope[i]*(FishLen-sims$L50[i])))

      # spawning duration with respect to fish length, for each set of simulated parameters
      SpawnDurEst[i,] = d50_2 - d50 # sim, fishlen

      for (k in 1:nFishLen) {
        P_t_s.sim[i,k,] = 1 / (1 + exp(-sims$slope1[i]*(DecDay-d50[k]))) *
          1 / (1 + exp(sims$slope2[i]*(DecDay-d50_2[k]))) * Height[k]

        # maximum proportion spawning with respect to fish length, for each set of simulated parameters
        MaxCurveHeightEst[i,k] = max(P_t_s.sim[i,k,]) # sim, fishlen, DecDay
      }

      cat("1:Resampling: i",i,"of",nsims,"nsims",'\n')
    } # nsims

    # relationship between maximum proportion spawning and fish size, with uncertainty
    MaxCurve.mean = apply(MaxCurveHeightEst[,], MARGIN=2, function(x) mean(x))
    MaxCurve.median = apply(MaxCurveHeightEst[,], MARGIN=2, function(x) quantile(x, 0.5))
    MaxCurve.lowCL = apply(MaxCurveHeightEst[,], MARGIN=2, function(x) quantile(x, 0.025))
    MaxCurve.uppCL = apply(MaxCurveHeightEst[,], MARGIN=2, function(x) quantile(x, 0.975))

    MaxPropSpawnVsFishLen = data.frame(FishLen=FishLen,
                                       MaxCurve.mean=MaxCurve.mean,
                                       MaxCurve.median=MaxCurve.median,
                                       MaxCurve.lowCL=MaxCurve.lowCL,
                                       MaxCurve.uppCL=MaxCurve.uppCL)

    # relationship between spawning duration and fish size, with uncertainty
    SpawnDur.mean = apply(SpawnDurEst[,], MARGIN=2, function(x) mean(x))
    SpawnDur.median = apply(SpawnDurEst[,], MARGIN=2, function(x) quantile(x, 0.5))
    SpawnDur.lowCL = apply(SpawnDurEst[,], MARGIN=2, function(x) quantile(x, 0.025))
    SpawnDur.uppCL = apply(SpawnDurEst[,], MARGIN=2, function(x) quantile(x, 0.975))
    SpawnDurVsFishLen = data.frame(FishLen=FishLen,
                                   SpawnDur.mean=SpawnDur.mean,
                                   SpawnDur.median=SpawnDur.median,
                                   SpawnDur.lowCL=SpawnDur.lowCL,
                                   SpawnDur.uppCL=SpawnDur.uppCL)
  } # if nsims > 0

  # get model curves with uncertainty for fish with respect to
  # midpoints of specified size categories
  lbnds = seq(MinLen,MaxLen-LenInc,LenInc)
  ubnds = lbnds + LenInc
  midpts = lbnds + (LenInc/2)
  nSizeClasses = length(midpts)
  DecDayPlot = seq(0,1,0.01)
  DecDay = DecDayPlot
  nDecDayPlot = length(DecDayPlot)
  ModelCurveEst <- array(dim=c(nsims,nSizeClasses,nDecDayPlot))
  if (nsims > 0) {
    for (i in 1:nsims) {
      for (k in 1:nSizeClasses) {
        fish_size= midpts[k]
        EstPar=as.numeric(sims[i,1:7])
        res = CalcDailySpawnProps_SpDurMod(fish_size, DecDay, EstPar)
        ModelCurveEst[i,k,] = res$P_t_s
      }
      cat("2:Resampling: i",i,"of",nsims,"nsims",'\n')
    }
  }

  # estimated curve describing proportion spawning, for fish of specified sizes,
  # corresponding to midpoints of specified length categories
  ModelCurvesVsFishLenCl <- array(dim=c(nSizeClasses,4,nDecDayPlot))
  dimnames(ModelCurvesVsFishLenCl)[[1]] <- midpts
  dimnames(ModelCurvesVsFishLenCl)[[2]] <- c("Mean","Median","Lw95CL","Up95CL")
  dimnames(ModelCurvesVsFishLenCl)[[3]] <- DecDay
  if (nsims > 0) {
    for (k in 1:nSizeClasses) {
      CurveEst.mean = apply(ModelCurveEst[,k,], MARGIN=2, function(x) mean(x))
      CurveEst.median = apply(ModelCurveEst[,k,], MARGIN=2, function(x) quantile(x, 0.5))
      CurveEst.lowCL = apply(ModelCurveEst[,k,], MARGIN=2, function(x) quantile(x, 0.025))
      CurveEst.uppCL = apply(ModelCurveEst[,k,], MARGIN=2, function(x) quantile(x, 0.975))
      ModelCurvesVsFishLenCl[k,1,]= CurveEst.mean
      ModelCurvesVsFishLenCl[k,2,]= CurveEst.median
      ModelCurvesVsFishLenCl[k,3,]= CurveEst.lowCL
      ModelCurvesVsFishLenCl[k,4,]= CurveEst.uppCL
      cat("3:Estimates from resampling: k",k,"of",nSizeClasses,"fish size classes",'\n')
    }
  }

  results = list(params=params,
             NLL=NLL,
             convergence=convergence,
             vcov.params=vcov.params,
             cor.params=cor.params,
             ParamEst=ParamEst,
             ParamEst.sim=ParamEst.sim,
             MaxPropSpawnVsFishLen=MaxPropSpawnVsFishLen,
             SpawnDurVsFishLen=SpawnDurVsFishLen,
             ModelCurvesVsFishLenCl=ModelCurvesVsFishLenCl,
             lbnds=lbnds,
             ubnds=ubnds,
             midpts=midpts,
             DecDayPlot=DecDayPlot)

  return(results)
}

#' Plot results from spawning duration model
#'
#' This function plots observed for the daily proportions of spawning capable fish for specified length categories,
#' and estimated daily proportions from spawning duration model (with associated 95 percent confidence limits),
#' (i.e. asymmetric double logistic model with variable height), and estimated relationships between proportion spawning
#' vs fish length at the estimated time of peak spawning, and of estimated spawning duration vs fish length
#'
#' @param DecDay decimal days of fish capture during the year
#' @param ObsSpawnDat data frame including observed fish lengths (FishLen), observed months of capture (MM)
#' and observed maturity status categories (ObsMatStatus, 0=immature, 1=mature)
#' @param MinLen Minimum fish length (for determining lengths to output expected model curves)
#' @param MaxLen maximum length of fish (for plotting data for length classes)
#' @param LenInc length increments of length classes, for plotting
#' @param lbnds lower bounds for length classes to be plotted
#' @param ubnds upper bounds for length classes to be plotted
#' @param params initial (transformed) values of model parameters, for the spawning duration model (including
#' lnL50, lnslope, logitPkSpawn, lnkappa, lnkappa2, lnslope1, lnslope2)
#' @param nsims number of resampled parameter values, used to calculate confidence intervals
#' @param FittedRes Outputs of fitted spawning duration model using the function GetSpawningDurationModelResults
#' @param PlotOpt 0=all plots, 1=prop spawning vs length, 2=spawn dur vs length, 3=prop spawn vs month (multiple plots, each by size category)
#'
#' @return plots associated with fit of spawning duration model
#' @examples
#' # Simulate spawning proportion data
#' set.seed(123)
#' GrowthEqn=1 # von Bertalanffy
#' nSamples = 5000
#' nSexes = 1
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = NA
#' Linf = 1000
#' vbK = 0.1
#' tzero = 0
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = 0.08
#' # Spawning duration model parameters
#' lnL50 = log(300)
#' lnslope= log(0.015)
#' logitPkSpawn = log(0.6/(1-0.6))
#' lnkappa = log(0.0002)
#' lnkappa2 = log(0.0002)
#' lnslope1 = log(30)
#' lnslope2 = log(30)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' Res=SimulateSpawningDurationData(params, GrowthEqn, nSamples, nSexes, MinAge, MaxAge,
#'                                  AgeStep, Ref_ages, Growth_params, Growth_cv)
#' ObsMatStatus=Res$ObsMatStatus
#' FishLen=Res$FishLen
#' DecDay=Res$DecDay
#' # Spawning duration model parameters - initial values
#' lnL50 = log(350)
#' lnslope= log(0.02)
#' logitPkSpawn = log(0.7/(1-0.7))
#' lnkappa = log(0.00015)
#' lnkappa2 = log(0.00035)
#' lnslope1 = log(25)
#' lnslope2 = log(35)
#' params = c(lnL50,lnslope,logitPkSpawn,lnkappa,lnkappa2,lnslope1,lnslope2)
#' nsims=20
#' MinLen=0
#' MaxLen=1100
#' LenInc=100
#' FittedRes=GetSpawningDurationModelResults(params, FishLen, DecDay, ObsMatStatus, MinLen, MaxLen, LenInc, nsims)
#' ObsSpawnDat = data.frame(FishLen=FishLen, ObsMatStatus=ObsMatStatus)
#' PlotOpt = 0 # 0=all plots, 1=prop spawning vs length, 2=spawn dur vs length, 3=prop spawn vs month (multiple plots, each by size category)
#' PlotSpawningDurationModelResults(DecDay, ObsSpawnDat, MinLen, MaxLen, LenInc, params, nsims, FittedRes, PlotOpt)
#' @export
PlotSpawningDurationModelResults <- function(DecDay, ObsSpawnDat, MinLen, MaxLen, LenInc, params, nsims, FittedRes, PlotOpt) {

  FishLen=ObsSpawnDat$FishLen
  ObsMatStatus=ObsSpawnDat$ObsMatStatus

  if (is.list(FittedRes)) {     # if model already fitted, can input results rather than refit
    Res =  FittedRes
  } else {
    Res=GetSpawningDurationModelResults(params, FishLen, DecDay, ObsMatStatus, MaxLen, LenInc, nsims)
  }

  if (PlotOpt==0) {
    par(mfrow=c(2,1), mar=c(4,4,2,2))
  } else {
    par(mfrow=c(1,1), mar=c(4,4,2,2))
  }

  # plot maximum proportion spawning vs fish length
  if (PlotOpt==0 | PlotOpt==1) {
    plot(Res$MaxPropSpawnVsFishLen[,1],Res$MaxPropSpawnVsFishLen[,3], "l", ylim=c(0,1), xlim=c(0,MaxLen), xaxt='n', yaxt="n",
         las=1, xlab=list("Fish length, mm",cex=1.2), ylab=list("Prop. Spawning",cex=1.2), cex.axis=1, bty='n')

    AddAxesAndTickLabelsToPlot(xmin=0, xmax=MaxLen, xint=LenInc, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=0,
                               lineval=0, lasval=2, xaxlabel = seq(0,MaxLen,LenInc), tcklen = 0.03)
    polygon(c(Res$MaxPropSpawnVsFishLen[,1],rev(Res$MaxPropSpawnVsFishLen[,1])),
            c(Res$MaxPropSpawnVsFishLen[,4],rev(Res$MaxPropSpawnVsFishLen[,5])),col="lightgrey",border="grey")
    lines(Res$MaxPropSpawnVsFishLen[,1],Res$MaxPropSpawnVsFishLen[,3])
  }

  # plot spawning duration vs fish length
  if (PlotOpt==0 | PlotOpt==2) {
    plot(Res$SpawnDurVsFishLen[,1],Res$SpawnDurVsFishLen[,3]*12, "l", ylim=c(0,12), xlim=c(0,MaxLen), xaxt='n', yaxt="n",
         las=1, xlab=list("Fish length, mm",cex=1.2), ylab=list("Spawning dur, months",cex=1.2), cex.axis=1, bty='n')
    AddAxesAndTickLabelsToPlot(xmin=0, xmax=MaxLen, xint=LenInc, ymin=0, ymax=12, yint=2, cexval=1,  cexaxisval=NA, lwdval=0,
                               lineval=0, lasval=2, xaxlabel = seq(0,MaxLen,LenInc), tcklen = 0.03)
    polygon(c(Res$SpawnDurVsFishLen[,1],rev(Res$SpawnDurVsFishLen[,1])),
            c(Res$SpawnDurVsFishLen[,4]*12,rev(Res$SpawnDurVsFishLen[,5]*12)),col="lightgrey",border="grey")
    lines(Res$SpawnDurVsFishLen[,1],Res$SpawnDurVsFishLen[,3]*12)
  }

  if (PlotOpt==0 | PlotOpt==3) {
    lbnds = seq(MinLen,MaxLen-LenInc,LenInc)
    ubnds = lbnds + LenInc
    midpts = lbnds + (LenInc/2)
    nLenCats = length(lbnds)
    floor(nLenCats/2)+1
    par(mfrow=c(3,2), mar=c(4,4,2,2))
    MMabb = substr(month.abb, 1, 1)
    for (i in 1:nLenCats) {
      subDat = ObsSpawnDat[ObsSpawnDat$FishLen >= lbnds[i] & ObsSpawnDat$FishLen < ubnds[i],]
      Probs = CalcMonthlyObsSpawnProps_SpDurMod(subDat)
      plot(seq(0.5,11.5,1), Probs, xaxt='n', yaxt="n", cex=0.8, cex.axis=1,bty='n', ylim=c(0,1), xlim=c(0,12),
           xlab = list(" Month",cex=1.2), ylab=list("Prop. spawning",cex=1.2), main=paste0(lbnds[i],"-",ubnds[i]," mm"), cex.main=1)
      AddAxesAndTickLabelsToPlot(xmin=0.5, xmax=11.5, xint=1, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=0,
                                 lineval=0, lasval=2, xaxlabel = MMabb[1:12], tcklen = 0.03)
      polygon(c(Res$DecDayPlot*12,rev(Res$DecDayPlot*12)),
              c(Res$ModelCurvesVsFishLen[i,3,],rev(Res$ModelCurvesVsFishLen[i,4,])),col="lightgrey",border="grey")
      lines(Res$DecDayPlot*12 ,Res$ModelCurvesVsFishLen[i,1,])
      points(seq(0.5,11.5,1), Probs)
    } # i
  }
}

#******************************
# Mixture distribution analyses
#******************************


#' Return expected proportions at length
#'
#' Return expected proportions at length, for mixture distribution analysis
#'
#' @keywords internal
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#' @return list object containing estimated proportion at age one, following feasible parameter value
#' check (EstPropOne), expected frequencies for first, second and third cohorts,
#' returning NA if not calculated due to specified number of cohorts (Expfreq_mode1, Expfreq_mode2, Expfreq_mode3),
#' EstPropZero = PropZero, expected overall proportion at length (ExpPropAtLen), expected overall
#' frequency at length (ExpFreqAtLen), values of penalties for normal disribution models, ensuring
#' mode 1 is not greater than mode 2, and mode 2 is not greater than mode 3
EstPropAtLen_MixtureDistn <- function(params) {

  if (length(which(c(2,4,5,6,8) == length(params))) == 0) {
    stop("Problem: not correct number of parameters")
  }

  # determine number of cohorts to be fitted, based on specified number of starting parameters
  if (length(params) == 2) {
    Cohorts=1
    if (DistnType == 1) { # normal
      Mean1 = exp(params[1])
      sd1 = exp(params[2])
    }
    if (DistnType == 2) { # gamma
      shape1 = exp(params[1])
      rate1 = exp(params[2])
    }
  }
  if (length(params) == 4) { # 2 cohorts with common sd or rate
    Cohorts=2
    if (DistnType == 1) { # normal
      Mean1 = exp(params[1])
      Mean2 = exp(params[2])
      sd1 = exp(params[3])
      sd2 = sd
    }
    if (DistnType == 2) { # gamma
      shape1 = exp(params[1])
      shape2 = exp(params[2])
      rate1 = exp(params[3])
      rate2 = rate1
    }
    PropZero = ilogit(params[4])
  }
  if (length(params) == 5) { # 2 cohorts with separate sds
    Cohorts=2
    if (DistnType == 1) { # normal
      Mean1 = exp(params[1])
      Mean2 = exp(params[2])
      sd1 = exp(params[3])
      sd2 = exp(params[4])
    }
    if (DistnType == 2) { # gamma
      shape1 = exp(params[1])
      shape2 = exp(params[2])
      rate1 = exp(params[3])
      rate2 = exp(params[4])
    }
    PropZero = ilogit(params[5])
  }
  if (length(params) == 6) { # 3 cohorts with common sds or rates
    Cohorts=3
    if (DistnType == 1) { # normal
      Mean1 = exp(params[1])
      Mean2 = exp(params[2])
      Mean3 = exp(params[3])
      sd1 = exp(params[4])
      sd2 = sd1
      sd3 = sd1
    }
    if (DistnType == 2) { # gamma
      shape1 = exp(params[1])
      shape2 = exp(params[2])
      shape3 = exp(params[3])
      rate1 = exp(params[4])
      rate2 = rate1
      rate3 = rate1
    }
    PropZero = ilogit(params[5])
    PropOne = ilogit(params[6])
  }
  if (length(params) == 8) { # 3 cohorts with separate sds or rates
    Cohorts=3
    if (DistnType == 1) { # normal
      Mean1 = exp(params[1])
      Mean2 = exp(params[2])
      Mean3 = exp(params[3])
      sd1 = exp(params[4])
      sd2 = exp(params[5])
      sd3 = exp(params[6])
    }
    if (DistnType == 2) { # gamma
      shape1 = exp(params[1])
      shape2 = exp(params[2])
      shape3 = exp(params[3])
      rate1 = exp(params[4])
      rate2 = exp(params[5])
      rate3 = exp(params[6])
    }
    PropZero = ilogit(params[7])
    PropOne = ilogit(params[8])
  }

  # get expected lengths for each cohort, and overall
  LbndSizeCl <- seq(MinSize, MaxSize-SizeInt, SizeInt)
  UbndSizeCl <- seq(SizeInt, MaxSize, SizeInt)
  Penalty1=0; Penalty2=0
  Expfreq_mode1=rep(NA,length(ObsFreq)); Expfreq_mode2=rep(NA,length(ObsFreq)); Expfreq_mode3=rep(NA,length(ObsFreq))

  if (Cohorts==1) {
    if (DistnType == 1) { # normal
      Expfreq_mode1 = (pnorm(UbndSizeCl, Mean1, sd1) -
                    pnorm(LbndSizeCl, Mean1, sd1)) * sum(ObsFreq)
    }
    if (DistnType == 2) { # gamma
      Expfreq_mode1 = (pgamma(UbndSizeCl, shape1, rate1) -
                    pgamma(LbndSizeCl, shape1, rate1)) * sum(ObsFreq)
    }
    sum_expfreq <- Expfreq_mode1 / sum(Expfreq_mode1)
  } # 1 cohort

  if (Cohorts==2) {

    # don't allow mean2 to be less than mean1 + 2mm
    if (DistnType == 1) { # normal
    if ((Mean2 - Mean1) < 2) {
      Penalty1 = 10 * (Mean2 - Mean1 - 2) ^ 2
      Mean2 = Mean1 + 2
    }
    }

    height1 <- sum(ObsFreq) * PropZero
    height2 <- sum(ObsFreq) * (1-PropZero)
    if (DistnType == 1) { # normal
      Expfreq_mode1 <- (pnorm(UbndSizeCl, Mean1, sd1) - pnorm(LbndSizeCl, Mean1, sd1)) * height1
      Expfreq_mode2 <- (pnorm(UbndSizeCl, Mean2, sd2) - pnorm(LbndSizeCl, Mean2, sd2)) * height2
    }
    if (DistnType == 2) { # gamma
      Expfreq_mode1 <- (pgamma(UbndSizeCl, shape1, rate1) - pgamma(LbndSizeCl, shape1, rate1)) * height1
      Expfreq_mode2 <- (pgamma(UbndSizeCl, shape2, rate2) - pgamma(LbndSizeCl, shape2, rate2)) * height2
    }
    sum_expfreq <- Expfreq_mode1 + Expfreq_mode2
  } # 2 cohorts

  if (Cohorts==3) {

    # don't allow mean2 to be less than mean1 + 2mm
    if (DistnType == 1) { # normal
    if ((Mean2 - Mean1) < 2) {
      Penalty1 = 10 * (Mean2 - Mean1 - 2) ^2
      Mean2 = Mean1 + 2
    }
    # don't allow mean3 to be less than mean2 + 2mm
    if ((Mean3 - Mean2) < 2) {
      Penalty1 = Penalty1 + (10 * (Mean3 - Mean2 - 2) ^2)
      Mean2 = Mean1 + 2
    }
    }
    # ensure sum of proportions is not > 1
    if ((PropZero + PropOne)  > 1.0) {
      Penalty2 = 10000 * ((PropZero + PropOne) - 1.0) ^ 2
      PropZero = PropZero / ((PropZero + PropOne))
      PropOne = 1 - PropZero
    }

    height1 <- sum(ObsFreq) * PropZero
    height2 <- sum(ObsFreq) * PropOne
    height3 <- sum(ObsFreq) * (1 - PropZero - PropOne)

    if (DistnType == 1) { # normal
      Expfreq_mode1 <- (pnorm(UbndSizeCl, Mean1, sd1) - pnorm(LbndSizeCl, Mean1, sd1)) * height1
      Expfreq_mode2 <- (pnorm(UbndSizeCl, Mean2, sd2) - pnorm(LbndSizeCl, Mean2, sd2)) * height2
      Expfreq_mode3 <- (pnorm(UbndSizeCl, Mean3, sd3) - pnorm(LbndSizeCl, Mean3, sd3)) * height3
    }
    if (DistnType == 2) { # gamma
      Expfreq_mode1 <- (pgamma(UbndSizeCl, shape1, rate1) - pgamma(LbndSizeCl, shape1, rate1)) * height1
      Expfreq_mode2 <- (pgamma(UbndSizeCl, shape2, rate2) - pgamma(LbndSizeCl, shape2, rate2)) * height2
      Expfreq_mode3 <- (pgamma(UbndSizeCl, shape3, rate3) - pgamma(LbndSizeCl, shape3, rate3)) * height3
    }
    sum_expfreq <- Expfreq_mode1 + Expfreq_mode2 + Expfreq_mode3
  } # 3 cohorts

  # get overall expected proportions and frequencies at length
  ExpPropAtLen = sum_expfreq / sum(sum_expfreq)
  ExpFreqAtLen = ExpPropAtLen * sum(ObsFreq)


  if (Cohorts==1) {
    results = list(EstPropZero = NA,
                   EstPropOne = NA,
                   Expfreq_mode1 = Expfreq_mode1,
                   Expfreq_mode2 = NA,
                   Expfreq_mode3 = NA,
                   ExpPropAtLen = ExpPropAtLen,
                   ExpFreqAtLen = ExpFreqAtLen,
                   Penalty1 = Penalty1,
                   Penalty2 = Penalty2)
  }
  if (Cohorts==2) {
    results = list(EstPropZero = PropZero,
                   EstPropOne = NA,
                   Expfreq_mode1 = Expfreq_mode1,
                   Expfreq_mode2 = Expfreq_mode2,
                   Expfreq_mode3 = NA,
                   ExpPropAtLen = ExpPropAtLen,
                   ExpFreqAtLen = ExpFreqAtLen,
                   Penalty1 = Penalty1,
                   Penalty2 = Penalty2)
  }
  if (Cohorts==3) {
    results = list(EstPropZero = PropZero,
                   EstPropOne = PropOne,
                   Expfreq_mode1 = Expfreq_mode1,
                   Expfreq_mode2 = Expfreq_mode2,
                   Expfreq_mode3 = Expfreq_mode3,
                   ExpPropAtLen = ExpPropAtLen,
                   ExpFreqAtLen = ExpFreqAtLen,
                   Penalty1 = Penalty1,
                   Penalty2 = Penalty2)
  }

  return(results)

} # end function


#' Calculate negative log-likelihood associated with the fit of a size mixture distribution model
#' and associated parameters values to size frequency data
#'
#' Calculates the negative log-likelihood associated with the fit of a size mixture
#' distribution model and associated parameters values to size frequency data,
#' with 1, 2 or 3 size modes. The model assumes each size mode is normally-distributed.
#' Four alternative models are available depending on the number of modes present in the
#' data, and assumptions regarding the spread of the various size modes.
#' Function requires size frequency data, and values for minimum size, maximum size
#' and size interval (stored in memory in R).
#' ObsFreq (Numeric Vector), MinSize, MaxSize, SizeInt (Numbers).
#' This function (with parameter inputs) can be passed into R optimisation routines (e.g. nlminb).
#'
#' @keywords internal
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#'
#' @return Negative log-likelihood associated with mixture model fit to size frequency data
CalcNLL_SizeMixtureDistn <- function(params) {

  # get expected proportions at length, for each cohort
  Res = EstPropAtLen_MixtureDistn(params)

  ExpPropAtLen = Res$ExpPropAtLen
  Penalty1 = Res$Penalty1
  Penalty2 = Res$Penalty2

  # calculate the multinomial likelihood
  NLL = -sum(ObsFreq * log(ExpPropAtLen + 0.0001)) + Penalty1 + Penalty2
  cat("NLL",NLL,"Penalty1",Penalty1,"Penalty2",Penalty2,'\n')

  results <- NLL

  return(results)

} # end function


#' Fit a size mixture distribution model to size frequency data
#'
#' This function fits a size mixture distribution model to size frequency data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' Function requires size frequency data (stored in memory in R), ObsFreq (Numeric Vector).
#'
#' @keywords internal
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#' @param DistnType 1=normal, 2=gamma
#' @param Cohorts 1, 2 or 3
#' @param MinSize minimum size
#' @param MaxSize maximum size
#' @param SizeInt size class interval
#' @param ObsFreq observed frequency data
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitMixtureDistnModel <- function(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq) {

  nlmb <- nlminb(params, CalcNLL_SizeMixtureDistn, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))
  results=nlmb
  return(results)
}

#' Get statistical outputs from a fitted size mixture distribution model
#'
#' This function fits a size mixture distribution model to size frequency data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' It provides various outputs in include convergence statistics, parameter estimated
#' and associated 95 calculated using the MASS package.
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#' @param DistnType 1=normal, 2=gamma
#' @param Cohorts 1, 2 or 3
#' @param MinSize minimum size
#' @param MaxSize maximum size
#' @param SizeInt size class interval
#' @param ObsFreq observed frequency data
#'
#' @return returns a list object containing negative log-likelihood (nll), nlminb convergence diagnostic
#' (convergence), sample size (SampleSize), estimates of model parameters with lower and upper 95
#' confidence limits (ParamEst), Estimates of standard deviations for estimated mean
#' (normal distribution assumption) or shape parameters (gamma distribution assumption) for
#' each fitted size mode (EstMean1_sd, EstMean2_sd, EstMean3_sd), point estimates
#' for model parameters (params) and variance-covariance matrix (vcov.params),
#' expected frequency for cohorts 1, 2 and 3 (Expfreq_mode1, Expfreq_mode2, Expfreq_mode3),
#' overall expected proportions at length across all cohorts (ExpPropAtLen), and
#' overall expected frequencies at length,across all cohorts (ExpFreqAtLen).
#' @examples
#' # Simulate data with 1 mode, specifying normal distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Cohorts=1
#' Mean1 = 20
#' sd1 = 5
#' ObsSize = round(rnorm(SampSize, Mean1, sd1),0)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' # Specify starting parameter values for Mean1 and sd1, and return associated negative log-likelihood
#' DistnType = 1 # 1=normal, 2=gamma
#' params = log(c(30, 5)) # normal
#' GetMixtureModelResults(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)
#' # Simulate data with 1 mode, specifying gamma distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Cohorts=1
#' Shape1 = 20
#' Rate1 = 1
#' ObsSize = round(rgamma(SampSize, Shape1, Rate1),0)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' DistnType = 2 # 1=normal, 2=gamma
#' params = log(c(20, 1)) # gamma
#' res=GetMixtureModelResults(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)
#' # Get mean and mode of fitted gamma distribution
#' alpha = res$ParamEst[1,1]
#' beta = 1 / res$ParamEst[2,1]
#' # Simulate data with 2 modes, specifying normal distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' # normal
#' Mean1 = 20
#' Mean2 = 35
#' sd1 = 5
#' sd2 = 5
#' PropZero = 0.7
#' ObsSize1 = round(rnorm(PropZero*SampSize, Mean1, sd1),0)
#' ObsSize2 = round(rnorm((1-PropZero)*SampSize, Mean2, sd2),0)
#' ObsSize = c(ObsSize1,ObsSize2)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' params = c(log(c(15, 40, 5, 5)), 0.5) # separate sd's for the 2 cohorts
#' Cohorts=2
#' DistnType = 1 # 1=normal, 2=gamma
#' GetMixtureModelResults(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)
#' # Simulate data with 2 modes, specifying gamma distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Shape1 = 20
#' Shape2 = 35
#' Rate1 = 1.5
#' Rate2 = 1
#' PropZero = 0.7
#' ObsSize1 = round(rgamma(PropZero*SampSize, Shape1, Rate1),0)
#' ObsSize2 = round(rgamma((1-PropZero)*SampSize, Shape2, Rate2),0)
#' ObsSize = c(ObsSize1,ObsSize2)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=TRUE)
#' ObsFreq = HistData$counts
#' params = c(log(c(15, 40, 1, 1)), 0.5) # separate sd's for the 2 cohorts
#' Cohorts=2
#' DistnType = 2 # 1=normal, 2=gamma
#' res=GetMixtureModelResults(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)
#' # Get means and modes of fitted gamma distributions
#' alpha = res$ParamEst[1,1]
#' beta = 1 / res$ParamEst[3,1]
#' Dist_mean1 = alpha * beta
#' Dist_mode1 = (alpha - 1) * beta
#' alpha = res$ParamEst[2,1]
#' beta = 1 / res$ParamEst[4,1]
#' Dist_mean2 = alpha * beta
#' Dist_mode2 = (alpha - 1) * beta
#' @export
GetMixtureModelResults <- function(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq) {

  if (length(which(c(2,4,5,6,8) == length(params))) == 0) {
    stop("Problem: not correct number of parameters")
  }

  # fit mixture model
  nlmb = FitMixtureDistnModel(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)

  # extract information useful for plotting
  params = nlmb$par
  Res=EstPropAtLen_MixtureDistn(params)
  Expfreq_mode1 = Res$Expfreq_mode1
  Expfreq_mode2 = Res$Expfreq_mode2
  Expfreq_mode3 = Res$Expfreq_mode3
  ExpPropAtLen = Res$ExpPropAtLen
  ExpFreqAtLen = Res$ExpFreqAtLen

  # calculate uncertainty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(nlmb$par, CalcNLL_SizeMixtureDistn)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  # back log-transform parameters
  if (length(params) == 2) { # 1 cohort
    if (DistnType == 1) { # normal
      EstMean1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstMean1_sd = ses[1]
      EstMean2_sd = NA
      EstMean3_sd = NA
      Estsd1 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      ParamEst = t(data.frame(EstMean1=round(EstMean1,2), Estsd1=round(Estsd1,2)))
    }
    if (DistnType == 2) { # gamma
      EstShape1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstShape1_sd = ses[1]
      EstShape2_sd = NA
      EstShape3_sd = NA
      EstRate1 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      ParamEst = t(data.frame(EstShape1=round(EstShape1,2), EstRate1=round(EstRate1,2)))
    }
  }
  if (length(params) == 4) { # 2 cohorts common sd
    if (DistnType == 1) { # normal
      EstMean1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstMean2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstMean1_sd = ses[1]
      EstMean2_sd = ses[2]
      EstMean3_sd = NA
      Estsd1 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      Estsd2 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstPropZero = ilogit(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      ParamEst = t(data.frame(EstMean1=round(EstMean1,2), EstMean2=round(EstMean2,2),
                              Estsd1=round(Estsd1,2), Estsd2=round(Estsd2,2),
                              EstPropZero=round(EstPropZero,2)))
      if (DistnType == 2) { # gamma
        EstShape1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
        EstShape2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstShape1_sd = ses[1]
        EstShape2_sd = ses[2]
        EstShape3_sd = NA
        EstRate1 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        EstRate2 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        EstPropZero = ilogit(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        ParamEst = t(data.frame(EstShape1=round(EstShape1,2), EstShape2=round(EstShape2,2),
                                EstRate1=round(EstRate1,2), Estsd2=round(EstRate2,2),
                                EstPropZero=round(EstPropZero,2)))
      }
    }
  }
  if (length(params) == 5) { # 2 cohorts, separate sds
    if (DistnType == 1) { # normal
      EstMean1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstMean2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstMean1_sd = ses[1]
      EstMean2_sd = ses[2]
      EstMean3_sd = NA
      Estsd1 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      Estsd2 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstPropZero = ilogit(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      ParamEst = t(data.frame(EstMean1=round(EstMean1,2), EstMean2=round(EstMean2,2),
                              Estsd1=round(Estsd1,2), Estsd2=round(Estsd2,2),
                              EstPropZero=round(EstPropZero,2)))
    }
    if (DistnType == 2) { # gamma
      EstShape1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstShape2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstShape1_sd = ses[1]
      EstShape2_sd = ses[2]
      EstShape3_sd = NA
      EstRate1 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstRate2 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstPropZero = ilogit(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      ParamEst = t(data.frame(EstShape1=round(EstShape1,2), EstShape2=round(EstShape2,2),
                              EstRate1=round(EstRate1,2), EstRate2=round(EstRate2,2),
                              EstPropZero=round(EstPropZero,2)))
    }
  }
  if (length(params) == 6) { # 3 cohorts common sd
    if (DistnType == 1) { # normal
      EstMean1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstMean2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstMean3 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstMean1_sd = ses[1]
      EstMean2_sd = ses[2]
      EstMean3_sd = ses[3]
      Estsd1 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      Estsd2 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      Estsd3 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstPropZero = ilogit(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      EstPropOne = ilogit(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
      ParamEst = t(data.frame(EstMean1=round(EstMean1,2), EstMean2=round(EstMean2,2),EstMean3=round(EstMean3,2),
                              Estsd1=round(Estsd1,2), Estsd2=round(Estsd2,2),Estsd3=round(Estsd3,2),
                              EstPropZero=round(EstPropZero,2),EstPropOne=round(EstPropOne,2)))
    }
    if (DistnType == 2) { # gamma
      EstShape1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstShape2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstShape3 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstShape1_sd = ses[1]
      EstShape2_sd = ses[2]
      EstShape3_sd = ses[3]
      EstRate1 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstRate2 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstRate3 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstPropZero = ilogit(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      EstPropOne = ilogit(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
      ParamEst = t(data.frame(EstShape1=round(EstShape1,2), EstShape2=round(EstShape2,2),EstShape3=round(EstShape3,2),
                              EstRate1=round(EstRate1,2), EstRate2=round(EstRate2,2),EstRate3=round(EstRate3,2),
                              EstPropZero=round(EstPropZero,2),EstPropOne=round(EstPropOne,2)))
    }
  }

  if (length(params) == 8) { # 3 cohorts
    if (DistnType == 1) { # normal
      EstMean1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstMean2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstMean3 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstMean1_sd = ses[1]
      EstMean2_sd = ses[2]
      EstMean3_sd = ses[3]
      Estsd1 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      Estsd2 = c(exp(nlmb$par[5]), exp(nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      Estsd3 = c(exp(nlmb$par[6]), exp(nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
      EstPropZero = ilogit(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7]))
      EstPropOne = ilogit(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8]))
      ParamEst = t(data.frame(EstMean1=round(EstMean1,2), EstMean2=round(EstMean2,2),EstMean3=round(EstMean3,2),
                              Estsd1=round(Estsd1,2), Estsd2=round(Estsd2,2),Estsd3=round(Estsd3,2),
                              EstPropZero=round(EstPropZero,2),EstPropOne=round(EstPropOne,2)))
    }
    if (DistnType == 2) { # gamma
      EstShape1 = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      EstShape2 = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      EstShape3 = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      EstShape1_sd = ses[1]
      EstShape2_sd = ses[2]
      EstShape3_sd = ses[3]
      EstRate1 = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      EstRate2 = c(exp(nlmb$par[5]), exp(nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
      EstRate3 = c(exp(nlmb$par[6]), exp(nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
      EstPropZero = ilogit(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7]))
      EstPropOne = ilogit(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8]))
      ParamEst = t(data.frame(EstShape1=round(EstShape1,2), EstShape2=round(EstShape2,2),EstShape3=round(EstShape3,2),
                              EstRate1=round(EstRate1,2), EstRate2=round(EstRate2,2),EstRate3=round(EstRate3,2),
                              EstPropZero=round(EstPropZero,2),EstPropOne=round(EstPropOne,2)))
    }
  }

  if (DistnType == 1) { # normal
    EstShape1_sd = NA
    EstShape2_sd = NA
    EstShape3_sd = NA
  }
  if (DistnType == 2) { # gamma
    EstMean1_sd = NA
    EstMean2_sd = NA
    EstMean3_sd = NA
  }

  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
  SampleSize = sum(ObsFreq)
  nll = nlmb$objective
  convergence = nlmb$convergence

  # get confidence limits for size frequencies
  Res=GetConfidenceLimitsForMixtureDistnCurve(params, vcov.params, DistnType)

  # store key results as summary
  ResultsSummary = list(ParamEst = ParamEst,
                        nll = nlmb$objective,
                        convergence = nlmb$convergence,
                        SampleSize = sum(ObsFreq))

  # store key model diagnostic information
  ModelDiag = list(EstMean1_sd = EstMean1_sd,
                   EstMean2_sd = EstMean2_sd,
                   EstMean3_sd = EstMean3_sd,
                   EstShape1_sd = EstShape1_sd,
                   EstShape2_sd = EstShape2_sd,
                   EstShape3_sd = EstShape3_sd,
                   gamma_mean1 = Res$gamma_mean1,
                   gamma_mean2 = Res$gamma_mean2,
                   gamma_mean3 = Res$gamma_mean3,
                   gamma_mode1 = Res$gamma_mode1,
                   gamma_mode2 = Res$gamma_mode2,
                   gamma_mode3 = Res$gamma_mode3,
                   gamma_mean1_lw = Res$gamma_mean1_lw,
                   gamma_mean2_lw = Res$gamma_mean2_lw,
                   gamma_mean3_lw = Res$gamma_mean3_lw,
                   gamma_mode1_lw = Res$gamma_mode1_lw,
                   gamma_mode2_lw = Res$gamma_mode2_lw,
                   gamma_mode3_lw = Res$gamma_mode3_lw,
                   gamma_mean1_up = Res$gamma_mean1_up,
                   gamma_mean2_up = Res$gamma_mean2_up,
                   gamma_mean3_up = Res$gamma_mean3_up,
                   gamma_mode1_up = Res$gamma_mode1_up,
                   gamma_mode2_up = Res$gamma_mode2_up,
                   gamma_mode3_up = Res$gamma_mode3_up)

  # store all results as a list object
  results = list(ParamEst = ParamEst,
                 params = nlmb$par,
                 nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 vcov.params = vcov.params,
                 cor.params = cor.params,
                 ses = ses,
                 Expfreq_mode1 = Expfreq_mode1,
                 Expfreq_mode2 = Expfreq_mode2,
                 Expfreq_mode3 = Expfreq_mode3,
                 ExpFreqAtLen = ExpFreqAtLen,
                 ExpPropAtLen = ExpPropAtLen,
                 ResultsSummary = ResultsSummary,
                 ModelDiag = ModelDiag)

  return(results)

}


#' Get expected size frequency from size mixture distribution model
#'
#' Gets expected size frequency from size mixture distribution model, used in calculating
#' confidence limits for mixture distribution curves applying a parametric resampling procedure
#'
#' @keywords internal
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#'
#' @return Expected size frequency at length (ExpFreqAtLen)
MixtureDistnCurve_SizeFreq <- function(params) {

  Res=EstPropAtLen_MixtureDistn(params)
  results = Res$ExpFreqAtLen

  return(results)
}


#' Get expected size frequency from size mixture distribution model
#'
#' Gets expected size frequency from size mixture distribution model, used in calculating
#' confidence limits for mixture distribution curves applying a parametric resampling procedure
#'
#' @keywords internal
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#' @param vcov.params estunated variance-covariance matrix for parameters of mixture model
#' @param DistnType 1=normal, 2=gamma
#' @return median values and lower and upper 95 percent confidence limits for expected size frequencies
#' at length calculated from parametric resampling procedure (sim.size.est, sim.size.low, sim.size.up),
#' median values and lower and upper 95 percent confidence limits for gamma mean and modal values,
#' from parametric resampling procedure, giving NAs if normal distribution specified, and if a parameter
#' is not required for specified model with only 1 or 2 modes (gamma_mean1, gamma_mean2,
#' gamma_mean3, gamma_mode1, gamma_mode2, gamma_mode3, gamma_mean1_lw, gamma_mean2_lw, gamma_mean3_lw,
#' gamma_mode1_lw, gamma_mode2_lw, gamma_mode2_lw, gamma_mean1_up, gamma_mean2_up, gamma_mean3_up,
#' gamma_mode1_up, gamma_mode2_up, gamma_mode3_up)
GetConfidenceLimitsForMixtureDistnCurve <- function(params, vcov.params, DistnType) {

  if (length(which(c(2,4,5,6,8) == length(params))) == 0) {
    stop("Problem: not correct number of parameters")
  }

  # store estimated parameter distributions
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = 500, params, vcov.params))
  # head(sims)

  # determine number of cohorts to be fitted, based on specified number of starting parameters
  if (length(params) == 2) { # Cohorts=1
    if (DistnType == 1) { # normal
      names(sims) = c("Mean1", "sd1")
    }
    if (DistnType == 2) { # gamma
      names(sims) = c("Shape1", "Rate1")
    }
  }
  if (length(params) == 4) { # Cohorts=2, common sd or rate
    if (DistnType == 1) { # normal
      names(sims) = c("Mean1", "Mean2","sd1","PropZero")
    }
    if (DistnType == 2) { # gamma
      names(sims) = c("Shape1", "Shape2","Rate1","PropZero")
    }
  }
  if (length(params) == 5) { # Cohorts=2, separate sds or rates
    if (DistnType == 1) { # normal
      names(sims) = c("Mean1", "Mean2","sd1","sd2","PropZero")
    }
    if (DistnType == 2) { # gamma
      names(sims) = c("Shape1", "Shape2","Rate1","Rate2","PropZero")
    }
  }
  if (length(params) == 6) { # Cohorts=3 common sd or rate
    if (DistnType == 1) { # normal
      names(sims) = c("Mean1", "Mean2","Mean3","sd1","PropZero","PropOne")
    }
    if (DistnType == 2) { # gamma
      names(sims) = c("Shape1", "Shape2","Shape3","Rate1","PropZero","PropOne")
    }
  }
  if (length(params) == 8) { # Cohorts=3 separate sds or rates
    if (DistnType == 1) { # normal
      names(sims) = c("Mean1", "Mean2","Mean3","sd1","sd2","sd3","PropZero","PropOne")
    }
    if (DistnType == 2) { # gamma
      names(sims) = c("Shape1", "Shape2","Shape3","Rate1","Rate2","Rate3","PropZero","PropOne")
    }
  }

  if (DistnType == 2) { # gamma
    gmean1 = rep(0,500)
    gmode1 = rep(0,500)
    gmean2 = rep(0,500)
    gmode2 = rep(0,500)
    gmean3 = rep(0,500)
    gmode3 = rep(0,500)
    for (i in 1:500) {
      if (length(params) == 2) { # Cohorts=1
        alpha1 = exp(sims[i,1])
        beta1 = 1 / exp(sims[i,2])
        gmean1[i] = alpha1 * beta1
        gmode1[i] = (alpha1 - 1) * beta1
      }

      if (length(params) == 4) { # Cohorts=2, common sd or rate
        alpha1 = exp(sims[i,1])
        alpha2 = exp(sims[i,2])
        beta1 = 1 / exp(sims[i,3])
        beta2 = beta1
        gmean1[i] = alpha1 * beta1
        gmode1[i] = (alpha1 - 1) * beta1
        gmean2[i] = alpha2 * beta2
        gmode2[i] = (alpha2 - 1) * beta2
      }

      if (length(params) == 5) { # Cohorts=2, separate sds or rates
        alpha1 = exp(sims[i,1])
        alpha2 = exp(sims[i,2])
        beta1 = 1 / exp(sims[i,3])
        beta2 = 1 / exp(sims[i,4])
        gmean1[i] = alpha1 * beta1
        gmode1[i] = (alpha1 - 1) * beta1
        gmean2[i] = alpha2 * beta2
        gmode2[i] = (alpha2 - 1) * beta2
      }
      if (length(params) == 6) { # Cohorts=3 common sd or rate
        alpha1 = exp(sims[i,1])
        alpha2 = exp(sims[i,2])
        alpha3 = exp(sims[i,3])
        beta1 = 1 / exp(sims[i,4])
        beta2 = beta1
        beta3 = beta1
        gmean1[i] = alpha1 * beta1
        gmode1[i] = (alpha1 - 1) * beta1
        gmean2[i] = alpha2 * beta2
        gmode2[i] = (alpha2 - 1) * beta2
        gmean3[i] = alpha3 * beta3
        gmode3[i] = (alpha3 - 1) * beta3
      }
      if (length(params) == 8) { # Cohorts=3 separate sds or rates
        alpha1 = exp(sims[i,1])
        alpha2 = exp(sims[i,2])
        alpha3 = exp(sims[i,3])
        beta1 = 1 / exp(sims[i,4])
        beta2 = 1 / exp(sims[i,5])
        beta3 = 1 / exp(sims[i,6])
        gmean1[i] = alpha1 * beta1
        gmode1[i] = (alpha1 - 1) * beta1
        gmean2[i] = alpha2 * beta2
        gmode2[i] = (alpha2 - 1) * beta2
        gmean3[i] = alpha3 * beta3
        gmode3[i] = (alpha3 - 1) * beta3
      }
    }
  }

  if (DistnType == 1) { # normal
    gamma_mean1 = NA
    gamma_mode1 = NA
    gamma_mean2 = NA
    gamma_mode2 = NA
    gamma_mean3 = NA
    gamma_mode3 = NA
    gamma_mean1_lw = NA
    gamma_mode1_lw = NA
    gamma_mean2_lw = NA
    gamma_mode2_lw = NA
    gamma_mean3_lw = NA
    gamma_mode3_lw = NA
    gamma_mean1_up = NA
    gamma_mode1_up = NA
    gamma_mean2_up = NA
    gamma_mode2_up = NA
    gamma_mean3_up = NA
    gamma_mode3_up = NA

  }
  if (DistnType == 2) { # gamma
    if (length(params) == 2) { # Cohorts=1
      gamma_mean1 = median(gmean1)
      gamma_mode1 = median(gmode1)
      gamma_mean2 = NA
      gamma_mode2 = NA
      gamma_mean3 = NA
      gamma_mode3 = NA
      gamma_mean1_lw = gamma_mean1 - (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_lw = gamma_mode1 - (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_lw = NA
      gamma_mode2_lw = NA
      gamma_mean3_lw = NA
      gamma_mode3_lw = NA
      gamma_mean1_up = gamma_mean1 + (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_up = gamma_mode1 + (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_up = NA
      gamma_mode2_up = NA
      gamma_mean3_up = NA
      gamma_mode3_up = NA
    }
    if (length(params) == 4 | length(params) == 5) { # Cohorts=2
      gamma_mean1 = median(gmean1)
      gamma_mode1 = median(gmode1)
      gamma_mean2 = median(gmean2)
      gamma_mode2 = median(gmode2)
      gamma_mean3 = NA
      gamma_mode3 = NA
      gamma_mean1_lw = gamma_mean1 - (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_lw = gamma_mode1 - (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_lw = gamma_mean2 - (1.96 * sqrt(sd(gmean2) / 500))
      gamma_mode2_lw = gamma_mode2 - (1.96 * sqrt(sd(gmode2) / 500))
      gamma_mean3_lw = NA
      gamma_mode3_lw = NA
      gamma_mean1_up = gamma_mean1 + (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_up = gamma_mode1 + (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_up = gamma_mean2 + (1.96 * sqrt(sd(gmean2) / 500))
      gamma_mode2_up = gamma_mode2 + (1.96 * sqrt(sd(gmode2) / 500))
      gamma_mean3_up = NA
      gamma_mode3_up = NA
    }
    if (length(params) == 6 | length(params) == 8) { # Cohorts=3
      gamma_mean1 = median(gmean1)
      gamma_mode1 = median(gmode1)
      gamma_mean2 = median(gmean2)
      gamma_mode2 = median(gmode2)
      gamma_mean3 = median(gmean3)
      gamma_mode3 = median(gmode3)
      gamma_mean1_lw = gamma_mean1 - (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_lw = gamma_mode1 - (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_lw = gamma_mean2 - (1.96 * sqrt(sd(gmean2) / 500))
      gamma_mode2_lw = gamma_mode2 - (1.96 * sqrt(sd(gmode2) / 500))
      gamma_mean3_lw = gamma_mean3 - (1.96 * sqrt(sd(gmean3) / 500))
      gamma_mode3_lw = gamma_mode3 - (1.96 * sqrt(sd(gmode3) / 500))
      gamma_mean1_up = gamma_mean1 + (1.96 * sqrt(sd(gmean1) / 500))
      gamma_mode1_up = gamma_mode1 + (1.96 * sqrt(sd(gmode1) / 500))
      gamma_mean2_up = gamma_mean2 + (1.96 * sqrt(sd(gmean2) / 500))
      gamma_mode2_up = gamma_mode2 + (1.96 * sqrt(sd(gmode2) / 500))
      gamma_mean3_up = gamma_mean3 + (1.96 * sqrt(sd(gmean3) / 500))
      gamma_mode3_up = gamma_mode3 + (1.96 * sqrt(sd(gmode3) / 500))
    }
  }

  sims.size = apply(X=sims[,], MARGIN=1, FUN=MixtureDistnCurve_SizeFreq)

  # Calculating the 2.5th an 97.5th percentile
  sim.size.est = apply(sims.size, 1, function(x) quantile(x, 0.5))
  sim.size.low = apply(sims.size, 1, function(x) quantile(x, 0.025))
  sim.size.up = apply(sims.size, 1, function(x) quantile(x, 0.975))

  results = list(sim.size.est = sim.size.est,
                 sim.size.low = sim.size.low,
                 sim.size.up = sim.size.up,
                 gamma_mean1 = gamma_mean1,
                 gamma_mean2 = gamma_mean2,
                 gamma_mean3 = gamma_mean3,
                 gamma_mode1 = gamma_mode1,
                 gamma_mode2 = gamma_mode2,
                 gamma_mode3 = gamma_mode3,
                 gamma_mean1_lw = gamma_mean1_lw,
                 gamma_mean2_lw = gamma_mean2_lw,
                 gamma_mean3_lw = gamma_mean3_lw,
                 gamma_mode1_lw = gamma_mode1_lw,
                 gamma_mode2_lw = gamma_mode2_lw,
                 gamma_mode3_lw = gamma_mode3_lw,
                 gamma_mean1_up = gamma_mean1_up,
                 gamma_mean2_up = gamma_mean2_up,
                 gamma_mean3_up = gamma_mean3_up,
                 gamma_mode1_up = gamma_mode1_up,
                 gamma_mode2_up = gamma_mode2_up,
                 gamma_mode3_up = gamma_mode3_up)

  return(results)

}


#' Plot fitted mixture model to size frequency data.
#'
#' @param params c(Mean1, sd1) normal distn: 1 mode, or
#' c(Shape1, Rate1) gamma distn: 1 mode, or
#' c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
#' c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
#' c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
#' c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
#' c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
#' c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds
#' @param DistnType 1=normal, 2=gamma
#' @param Cohorts 1, 2 or 3
#' @param ObsFreq observed frequency data
#' @param MinSize minimum size
#' @param MaxSize maximum size
#' @param SizeInt size class interval
#' @param ymax y axis maximum
#' @param xmax x axis maximum
#' @param yint interval for y axis
#' @param xint interval for x axis
#' @param GraphTitle title for graph
#' @param xaxis_lab label for x axis
#' @param yaxis_lab label for y axis
#' @param set.plot.par logical, true sets program default
#' @param PlotHist logical, true is histogram, false is scatter plot
#' @param PlotCLs logical, true plots 95 percent confidence limits
#'
#' @return fitted mixture distribution curves either scatter plot or histogram
#'
#' @examples
#' # Simulate data with 1 mode, specifying normal distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Cohorts=1
#' Mean1 = 20
#' sd1 = 5
#' ObsSize = round(rnorm(SampSize, Mean1, sd1),0)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' # Specify starting parameter values for Mean1 and sd1, and return associated negative log-likelihood
#' DistnType = 1 # 1=normal, 2=gamma
#' params = log(c(30, 5)) # normal
#' PlotMixtureDistnResults(params, DistnType, Cohorts, ObsFreq, MinSize, MaxSize, SizeInt, ymax=100, xmax=50, yint=20, xint=10,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, set.plot.par=TRUE, PlotHist=TRUE, PlotCLs=FALSE)
#' # Simulate data with 1 mode, specifying gamma distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Cohorts=1
#' Shape1 = 20
#' Rate1 = 1
#' ObsSize = round(rgamma(SampSize, Shape1, Rate1),0)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' DistnType = 2 # 1=normal, 2=gamma
#' params = log(c(20, 1)) # gamma
#' PlotMixtureDistnResults(params, DistnType, Cohorts, ObsFreq, MinSize, MaxSize, SizeInt, ymax=100, xmax=50, yint=20, xint=10,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, set.plot.par=TRUE, PlotHist=TRUE, PlotCLs=FALSE)
#' # Simulate data with 2 modes, specifying normal distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' # normal
#' Mean1 = 20
#' Mean2 = 35
#' sd1 = 5
#' sd2 = 5
#' PropZero = 0.7
#' ObsSize1 = round(rnorm(PropZero*SampSize, Mean1, sd1),0)
#' ObsSize2 = round(rnorm((1-PropZero)*SampSize, Mean2, sd2),0)
#' ObsSize = c(ObsSize1,ObsSize2)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=FALSE)
#' ObsFreq = HistData$counts
#' params = c(log(c(15, 40, 5, 5)), 0.5) # separate sd's for the 2 cohorts
#' Cohorts=2
#' DistnType = 1 # 1=normal, 2=gamma
#' PlotMixtureDistnResults(params, DistnType, Cohorts, ObsFreq, MinSize, MaxSize, SizeInt, ymax=100, xmax=50, yint=20, xint=10,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, set.plot.par=TRUE, PlotHist=TRUE, PlotCLs=FALSE)
#' # Simulate data with 2 modes, specifying gamma distribution
#' set.seed(123)
#' MinSize = 0
#' MaxSize = 60
#' SizeInt = 1
#' SampSize = 1000
#' Shape1 = 20
#' Shape2 = 35
#' Rate1 = 1.5
#' Rate2 = 1
#' PropZero = 0.7
#' ObsSize1 = round(rgamma(PropZero*SampSize, Shape1, Rate1),0)
#' ObsSize2 = round(rgamma((1-PropZero)*SampSize, Shape2, Rate2),0)
#' ObsSize = c(ObsSize1,ObsSize2)
#' HistData=hist(ObsSize, breaks=seq(0,MaxSize,SizeInt), right=FALSE, plot=TRUE)
#' ObsFreq = HistData$counts
#' params = c(log(c(15, 40, 1, 1)), 0.5) # separate sd's for the 2 cohorts
#' Cohorts=2
#' DistnType = 2 # 1=normal, 2=gamma
#' PlotMixtureDistnResults(params, DistnType, Cohorts, ObsFreq, MinSize, MaxSize, SizeInt, ymax=100, xmax=50, yint=20, xint=10,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, set.plot.par=TRUE, PlotHist=TRUE, PlotCLs=FALSE)
#' @export
PlotMixtureDistnResults <- function(params, DistnType, Cohorts, ObsFreq, MinSize, MaxSize, SizeInt, ymax, xmax, yint, xint,
                                    GraphTitle, xaxis_lab, yaxis_lab, set.plot.par, PlotHist, PlotCLs) {

  LbndSizeCl <- seq(0, MaxSize-SizeInt, SizeInt)
  MidPtSizeCl <- LbndSizeCl + (SizeInt/2)
  SizeIndivs=rep(LbndSizeCl,ObsFreq)
  if (xmax<max(SizeIndivs)) xmax = trunc(max(SizeIndivs)/xint)*xint+xint

  # get default axis limits and intervals
  xlims=Get_xaxis_scale(0:MaxSize)
  ylims=Get_yaxis_scale(ObsFreq)

  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "Size, mm"
  if (is.na(yaxis_lab)) yaxis_lab = "Frequency"

  .pardefault <- par(no.readonly = TRUE) # store default par settings

  if (set.plot.par == T) {
    par(mfrow=c(1,1), mar=c(5,4,2,2), oma=c(2,2,2,2))
  }

  # get parameter estimates
  Res=FitMixtureDistnModel(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)
  params = Res$par
  Res=EstPropAtLen_MixtureDistn(params)

  if (PlotHist == FALSE) {
    plot(MidPtSizeCl, ObsFreq, "p", xlim=c(0,xmax), ylim=c(0,ymax), cex=0.5, pch=16,
         frame=F, xaxt = 'n', yaxt = 'n', xlab="", ylab="", cex.main=1, main=GraphTitle)
  }
  if (PlotHist == TRUE) {
    HistData=hist(SizeIndivs, breaks=0:xmax, right=FALSE, plot=TRUE, col="light grey", border="black",
                  main=GraphTitle, cex.main=1, xaxt = 'n', yaxt = 'n', xlab="", ylab="", ylim=c(0,ymax))
  }

  lines(MidPtSizeCl, Res$ExpFreqAtLen, lwd=2, lty="solid") # overall model
  if (!is.na(sum(Res$Expfreq_mode1))) {
    lines(MidPtSizeCl, Res$Expfreq_mode1, lwd=1, lty="dotted") # 1st cohort
  }
  if (!is.na(sum(Res$Expfreq_mode2))) {
    lines(MidPtSizeCl, Res$Expfreq_mode2, lwd=1, lty="dotted") # 2nd cohort
  }
  if (!is.na(sum(Res$Expfreq_mode3))) {
    lines(MidPtSizeCl, Res$Expfreq_mode3, lwd=1, lty="dotted") # 3rd cohort
  }
  AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=NA,  cexaxisval=1, lwdval=0,
                             lineval=-0.2, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
  mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2)
  mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2)
  if (length(params) == 2) Modes=1
  if (length(params) == 4) Modes=2 # common sd
  if (length(params) == 5) Modes=2 # separate sds
  if (length(params) == 6) Modes=3 # common sd
  if (length(params) == 8) Modes=3 # separate sds
  mtext(paste("Modes = ",Modes, sep=""), adj=0.1,side=3,las=1,line=-1,cex=0.8)

  if (PlotCLs == TRUE) {

    # fit mixture model and get results, including variance-covariance matrix
    Res=GetMixtureModelResults(params, DistnType, Cohorts, MinSize, MaxSize, SizeInt, ObsFreq)

    # get get parameter estimates and variance-covariance matrix
    params=Res$params
    vcov.params=Res$vcov.params

    # get confidence limits for size frequencies
    Res=GetConfidenceLimitsForMixtureDistnCurve(params, vcov.params, DistnType)

    # plot confidence limits
    if (PlotHist == TRUE) {
      x = c(MidPtSizeCl,rev(MidPtSizeCl)) # using shading for 95% CLs
      y = c(Res$sim.size.low, rev(Res$sim.size.up))
      polygon(x,y,col="light grey",border=NA)
      lines(MidPtSizeCl, Res$sim.size.est, "l", lty="solid")
      lines(MidPtSizeCl, Res$sim.size.low, "l", lty="dotted")
      lines(MidPtSizeCl, Res$sim.size.up, "l", lty="dotted")
    }
    if (PlotHist == FALSE) {
      x = c(MidPtSizeCl,rev(MidPtSizeCl)) # using shading for 95% CLs
      y = c(Res$sim.size.low, rev(Res$sim.size.up))
      polygon(x,y,col="light grey",border=NA)
      lines(MidPtSizeCl, Res$sim.size.est, "l", lty="solid")
      lines(MidPtSizeCl, Res$sim.size.low, "l", lty="dotted")
      lines(MidPtSizeCl, Res$sim.size.up, "l", lty="dotted")
      points(MidPtSizeCl, ObsFreq,cex=0.5, pch=16)
    }
  }

  # reset default par options
  par(.pardefault)
}


#****************************
# GROWTH - length-at-age data
#****************************



#' Simulate fish length at age data from von Bertalanffy growth curve
#'
#' Simulates fish length at age data for individual fish, and mean length at age data
#'
#' @param GrowthEqn 1=von Bertalanffy, 2=Schnute, 3=Somers seasonal curve, 4=von Bertalanffy with sex or area-divergent growth
#' @param nSamples number of samples
#' @param nSexes number of sexes 1=single or combined sex, 2=separate sexes
#' @param MinAge  minimum age
#' @param MaxAge maximum age
#' @param AgeStep age interval
#' @param Linf asymptotic length
#' @param vbK von Bertalanffy growth coefficient
#' @param tzero hypothetical length at age zero
#' @param Growth_sd standard deviation of lengths at age
#' @return lengths for individual fish (ObsLen), age classes in sample (ObsAgeCl),
#' mean length for each age class for females and males or combined sexes (ObsMeanLen,
#' FemObsMeanLen, MalObsMeanLen), standard deviations for mean
#' lengths at age for females and males or combined sexes (ObsMeanLensd, FemObsMeanLensd,
#' MalObsMeanLensd), lower and upper 95 percent confidence limits for
#' mean lengths at age for females and males or combined sexes (ObsMeanLen_lw95CL, ObsMeanLen_hi95CL,
#' FemObsMeanLen_lw95CL, MalObsMeanLen_lw95CL, FemObsMeanLen_hi95CL, MalObsMeanLen_hi95CL)
#'
#' @examples
#' set.seed(123)
#' # 1 sex
#' GrowthEqn=1
#' nSamples = 1000
#' nSexes = 1
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = NA
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = 0.08
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # plot lengths at age
#' plot(Res$ObsAge, Res$ObsLen, ylim=c(0,600))
#' # plot mean length at age data
#' ObsAgeCl=Res$ObsAgeCl
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLen_lw95CL=Res$ObsMeanLen_lw95CL
#' ObsMeanLen_hi95CL=Res$ObsMeanLen_hi95CL
#' plot(ObsAgeCl, ObsMeanLen, ylim=c(0,500))
#' arrows(ObsAgeCl, ObsMeanLen_lw95CL, ObsAgeCl, ObsMeanLen_hi95CL,
#'        code=3, angle=90,length=0.02, col='black')
#' # # 2 sexes
#' # GrowthEqn = 1
#' # nSamples = c(300,300)
#' # nSexes = 2
#' # MinAge = 1
#' # MaxAge = 20
#' # AgeStep = 1
#' # Linf = c(300,250)
#' # vbK = c(0.3,0.3)
#' # tzero = c(0,0)
#' # Growth_params = c(Linf,vbK,tzero)
#' # Growth_cv = c(0.08,0.08)
#' # Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # # plot lengths at age
#' # for (i in 1:nSexes) {
#' #   ObsAge=as.vector(unlist(Res$ObsAge[i,]))
#' #   ObsLen=as.vector(unlist(Res$ObsLen[i,]))
#' #   if (i==1) {
#' #     plot(ObsAge, ObsLen, ylim=c(0,500))
#' #   } else {
#' #     points(ObsAge, ObsLen, col="blue")
#' #   }
#' # }
#' # # plot mean lengths at age
#' # ObsAgeCl=Res$ObsAgeCl
#' # ObsMeanLen=Res$FemObsMeanLen
#' # ObsMeanLen_lw95CL=Res$FemObsMeanLen_lw95CL
#' # ObsMeanLen_hi95CL=Res$FemObsMeanLen_hi95CL
#' # plot(ObsAgeCl, ObsMeanLen, ylim=c(0,500))
#' # arrows(ObsAgeCl, ObsMeanLen_lw95CL, ObsAgeCl, ObsMeanLen_hi95CL,
#' #        code=3, angle=90,length=0.02, col='black')
#' # ObsMeanLen=Res$MalObsMeanLen
#' # ObsMeanLen_lw95CL=Res$MalObsMeanLen_lw95CL
#' # ObsMeanLen_hi95CL=Res$MalObsMeanLen_hi95CL
#' # points(ObsAgeCl, ObsMeanLen, col="blue")
#' # arrows(ObsAgeCl, ObsMeanLen_lw95CL, ObsAgeCl, ObsMeanLen_hi95CL,
#' #        code=3, angle=90,length=0.02, col='blue')
#' # # von Bertalanffy growth equation with divergent growth
#' # GrowthEqn = 4
#' # nSexes = 2 # single or combined sex
#' # nSamples = 500
#' # MinAge = 1
#' # MaxAge = 20
#' # AgeStep = 1
#' # Linf = c(300,500)
#' # vbK = c(0.3,0.3)
#' # tzero1 = 0
#' # tdiverge = 3
#' # Growth_params = c(Linf, vbK, tzero1, tdiverge)
#' # Ref_ages = NA
#' # Growth_cv = 0.08
#' # Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # x=which(Res$ObsDat$ObsSex==1)
#' # plot(Res$ObsAge[x], Res$ObsLen[x])
#' # y=which(Res$ObsDat$ObsSex==2)
#' # points(Res$ObsAge[y], Res$ObsLen[y],col="blue")
#' @export
SimulateLengthAtAgeData <- function(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep,
                                    Ref_ages, Growth_params, Growth_cv) {

  if (nSexes==1 & GrowthEqn<4) {
    ObsAge = runif(nSamples, MinAge, MaxAge)

    if (GrowthEqn == 1) { # von Bertalanffy
      Linf=Growth_params[1]
      vbK=Growth_params[2]
      tzero=Growth_params[3]
      MeanLen = Linf * (1 - exp(-vbK * (ObsAge - tzero)))
    }
    if (GrowthEqn == 2) { # Schnute
      t1=Ref_ages[1]; t2=Ref_ages[2]
      y1=Growth_params[1]; y2=Growth_params[2]
      a=Growth_params[3]; b=Growth_params[4]
      MeanLen <- rep(NA,length(ObsAge))
      for (i in 1:nSamples) {
        Age <- ObsAge[i]
        MeanLen[i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
      }
    }
    if (GrowthEqn == 3) { # Somers seasonal curve
      Linf = Growth_params[1]
      vbK = Growth_params[2]
      tzero = Growth_params[3]
      tc = Growth_params[4]
      C = Growth_params[5]
      S_age_t = sin(2 * pi * (ObsAge - tc))
      S_tzero = sin(2 * pi * (tzero - tc))
      MeanLen = Linf * (1 - exp(-vbK * (ObsAge - tzero +
                                          (C / (2 * pi)) * (S_age_t - S_tzero))))
    }

    # get random lengths at age
    ObsLen = rnorm(nSamples, MeanLen, Growth_cv * MeanLen)

    # observed mean length and sd data
    ObsDat = data.frame(ObsAge=ObsAge, ObsLen=ObsLen)
    ObsDat$AgeCl = trunc(ObsAge/AgeStep)*AgeStep
    MeanLenAtAge = aggregate(ObsLen ~ AgeCl, data=ObsDat, mean)
    ObsAgeCl = MeanLenAtAge[,1]
    ObsMeanLen = MeanLenAtAge[,2]
    Dat2 = aggregate(ObsLen ~ AgeCl, data=ObsDat, sd)
    Dat3 = aggregate(ObsLen ~ AgeCl, data=ObsDat, length)
    ObsMeanLensd = Dat2[,2]; ObsMeanLennObs = Dat3[,2]
    ObsMeanLense = ObsMeanLensd / sqrt(ObsMeanLennObs-1)
    ObsMeanLen_lw95CL = ObsMeanLen - (1.96 * ObsMeanLense); ObsMeanLen_hi95CL = ObsMeanLen + (1.96 * ObsMeanLense)
    ObsMeanLen_lw95PL = ObsMeanLen - (1.96 * ObsMeanLensd); ObsMeanLen_hi95PL = ObsMeanLen + (1.96 * ObsMeanLensd)
    FemObsMeanLen = NA; FemObsMeanLensd = NA; FemObsMeanLennObs = NA; FemObsMeanLense = NA
    FemObsMeanLen_lw95CL = NA; FemObsMeanLen_hi95CL = NA; FemObsMeanLen_lw95PL = NA; FemObsMeanLen_hi95PL = NA
    MalObsMeanLen = NA; MalObsMeanLensd = NA; MalObsMeanLennObs = NA; MalObsMeanLense = NA
    MalObsMeanLen_lw95CL = NA; MalObsMeanLen_hi95CL = NA; MalObsMeanLen_lw95PL = NA; MalObsMeanLen_hi95PL = NA
  }

  if (nSexes==2 & GrowthEqn<4) { # 2 sexes (or areas of same sex)
    ObsAge = data.frame(matrix(nrow=2,ncol=max(nSamples)))
    colnames(ObsAge) = 1:max(nSamples)
    MeanLen = ObsAge; ObsLen = ObsAge

    for (i in 1:nSexes) {
      tempObsAge = runif(nSamples[i], MinAge, MaxAge)
      # get mean lengths
      if (GrowthEqn == 1) { # von Bertalanffy
        Linf=Growth_params[1:2]
        vbK=Growth_params[3:4]
        tzero=Growth_params[5:6]
        tempMeanLen = Linf[i] * (1 - exp(-vbK[i] * (tempObsAge - tzero[i])))
      }
      if (GrowthEqn == 2) { # Schnute
        t1=Ref_ages[1]; t2=Ref_ages[2]
        if (i==1) {
          y1=Growth_params[1]; y2=Growth_params[3]
          a=Growth_params[5]; b=Growth_params[7]
        } else {
          y1=Growth_params[2]; y2=Growth_params[4]
          a=Growth_params[6]; b=Growth_params[8]
        }
        tempMeanLen = rep(0,length(nSamples[i]))
        for (j in 1:nSamples[i]) {
          Age <- tempObsAge[j]
          tempMeanLen[j] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
        }
      }
      if (GrowthEqn == 3) { # Somers seasonal curve
        Linf = Growth_params[1:2]
        vbK = Growth_params[3:4]
        tzero = Growth_params[5:6]
        tc = Growth_params[7:8]
        C = Growth_params[9:10]
        S_age_t = sin(2 * pi * (tempObsAge - tc[i]))
        S_tzero = sin(2 * pi * (tzero[i] - tc[i]))
        tempMeanLen = Linf[i] * (1 - exp(-vbK[i] * (tempObsAge - tzero[i] +
                                                      (C[i] / (2 * pi)) * (S_age_t - S_tzero))))
      }

      # get random lengths at age
      tempObsLen = rnorm(nSamples[i], tempMeanLen, Growth_cv[i] * tempMeanLen)

      # observed mean length, sd and se, 95 percent confidence limits and prediction limits
      ObsDat = data.frame(tempObsAge=tempObsAge,
                          tempObsLen=tempObsLen)
      ObsDat$tempAgeCl = trunc(tempObsAge/AgeStep)*AgeStep
      MeanLenAtAge = aggregate(tempObsLen ~ tempAgeCl, data=ObsDat, mean)
      ObsAgeCl = MeanLenAtAge[,1]
      ObsMeanLen = MeanLenAtAge[,2]
      Dat2 = aggregate(tempObsLen ~ tempAgeCl, data=ObsDat, sd)
      Dat3 = aggregate(tempObsLen ~ tempAgeCl, data=ObsDat, length)
      ObsMeanLensd = Dat2[,2]; ObsMeanLennObs = Dat3[,2]
      ObsMeanLense = ObsMeanLensd / sqrt(ObsMeanLennObs-1)
      ObsMeanLen_lw95CL = ObsMeanLen - (1.96 * ObsMeanLense); ObsMeanLen_hi95CL = ObsMeanLen + (1.96 * ObsMeanLense)
      ObsMeanLen_lw95PL = ObsMeanLen - (1.96 * ObsMeanLensd); ObsMeanLen_hi95PL = ObsMeanLen + (1.96 * ObsMeanLensd)

      if (i==1) {
        FemObsMeanLen=MeanLenAtAge[,2]; FemObsMeanLensd=ObsMeanLensd
        FemObsMeanLennObs=ObsMeanLennObs; FemObsMeanLense=ObsMeanLense
        FemObsMeanLen_lw95CL=ObsMeanLen_lw95CL; FemObsMeanLen_hi95CL=ObsMeanLen_hi95CL
        FemObsMeanLen_lw95PL=ObsMeanLen_lw95PL; FemObsMeanLen_hi95PL=ObsMeanLen_hi95PL
      } else {
        MalObsMeanLen = MeanLenAtAge[,2]; MalObsMeanLensd = ObsMeanLensd
        MalObsMeanLennObs=ObsMeanLennObs; MalObsMeanLense=ObsMeanLense
        MalObsMeanLen_lw95CL=ObsMeanLen_lw95CL; MalObsMeanLen_hi95CL=ObsMeanLen_hi95CL
        MalObsMeanLen_lw95PL=ObsMeanLen_lw95PL; MalObsMeanLen_hi95PL=ObsMeanLen_hi95PL
      }
      ObsAge[i,1:nSamples[i]] = tempObsAge
      ObsLen[i,1:nSamples[i]] = tempObsLen
    } # sex
    ObsMeanLen=NA; ObsMeanLensd=NA; ObsMeanLennObs=NA; ObsMeanLense=NA
    ObsMeanLen_lw95CL=NA; ObsMeanLen_hi95CL=NA; ObsMeanLen_lw95PL=NA; ObsMeanLen_hi95PL=NA
  } # i

  if (GrowthEqn == 4) { # von Bertalanffy with divergent growth

    MeanLen <- rep(NA,length(nSamples))
    ObsAge = runif(nSamples, MinAge, MaxAge)
    ObsSex = sample(c(1,2), nSamples, replace = TRUE)

    Linf = Growth_params[1:2]
    if (length(Growth_params)==6) { # separate vbk
      vbK = Growth_params[3:4]
      tzero1 = Growth_params[5]
      tdiverge = Growth_params[6]
    }

    if (length(Growth_params)==5) { # common vbK
      vbK = rep(0,2)
      vbK[1] = Growth_params[3]; vbK[2] = Growth_params[3]
      tzero1 = Growth_params[4]
      tdiverge = Growth_params[5]
    }

    Ldiverge = Linf[1] * (1 - exp(-vbK[1] * (tdiverge - tzero1)))
    if (Ldiverge > Linf[2] - 1)  Ldiverge = Linf[2] - 1
    tzero2 = tzero1 - (1/vbK[1])*log(1-(Ldiverge/Linf[1])) +
      (1/vbK[2])*log(1-(Ldiverge/Linf[2]))

    for (j in 1:nSamples) {
      sex = ObsSex[j]
      if (sex == 1 | ObsAge[j] <= tdiverge) {
        MeanLen[j] = Linf[1] * (1 - exp(-vbK[1] * (ObsAge[j] - tzero1)))
      }
      if (sex == 2 & ObsAge[j] > tdiverge)
        MeanLen[j] = Linf[2] * (1 - exp(-vbK[2] * (ObsAge[j] - tzero2)))
    }

    # get random lengths at age
    ObsLen = rnorm(nSamples, MeanLen, Growth_cv * MeanLen)

    # observed mean length, sd and se, 95 percent confidence limits and prediction limits
    ObsDat = data.frame(ObsAge=ObsAge, ObsLen=ObsLen, ObsSex=ObsSex)
    ObsDat$AgeCl = trunc(ObsAge/AgeStep)*AgeStep
    MeanLenAtAge = aggregate(ObsLen ~ AgeCl, data=ObsDat, mean)
    ObsAgeCl = MeanLenAtAge[,1]
    ObsMeanLen = MeanLenAtAge[,2]
    ObsMeanLen=NA; ObsMeanLensd=NA; ObsMeanLennObs=NA; ObsMeanLennObs=NA
    ObsMeanLen_lw95CL=NA; ObsMeanLen_hi95CL=NA; ObsMeanLen_lw95PL=NA; ObsMeanLen_hi95PL=NA
    FemObsMeanLen=NA; FemObsMeanLensd=NA; FemObsMeanLennObs=NA; FemObsMeanLense=NA
    FemObsMeanLen_lw95CL=NA; FemObsMeanLen_hi95CL=NA; FemObsMeanLen_lw95PL=NA; FemObsMeanLen_hi95PL=NA
    MalObsMeanLen=NA; MalObsMeanLensd=NA; MalObsMeanLensd=NA; MalObsMeanLennObs=NA
    MalObsMeanLen_lw95CL=NA; MalObsMeanLen_hi95CL=NA; MalObsMeanLen_lw95PL=NA; MalObsMeanLen_hi95PL=NA
  }

  results = list(ObsDat=ObsDat,
                 ObsAge=ObsAge,
                 ObsLen=ObsLen,
                 ObsAgeCl=ObsAgeCl,
                 ObsMeanLen=ObsMeanLen,
                 ObsMeanLensd=ObsMeanLensd,
                 ObsMeanLennObs=ObsMeanLennObs,
                 ObsMeanLense=ObsMeanLense,
                 ObsMeanLen_lw95CL=ObsMeanLen_lw95CL,
                 ObsMeanLen_hi95CL=ObsMeanLen_hi95CL,
                 ObsMeanLen_lw95PL=ObsMeanLen_lw95PL,
                 ObsMeanLen_hi95PL=ObsMeanLen_hi95PL,
                 FemObsMeanLen=FemObsMeanLen,
                 FemObsMeanLensd=FemObsMeanLensd,
                 FemObsMeanLennObs=FemObsMeanLennObs,
                 FemObsMeanLense=FemObsMeanLense,
                 FemObsMeanLen_lw95CL=FemObsMeanLen_lw95CL,
                 FemObsMeanLen_hi95CL=FemObsMeanLen_hi95CL,
                 FemObsMeanLen_lw95PL=FemObsMeanLen_lw95PL,
                 FemObsMeanLen_hi95PL=FemObsMeanLen_hi95PL,
                 MalObsMeanLen=MalObsMeanLen,
                 MalObsMeanLensd=MalObsMeanLensd,
                 MalObsMeanLennObs=MalObsMeanLennObs,
                 MalObsMeanLense=MalObsMeanLense,
                 MalObsMeanLen_lw95CL=MalObsMeanLen_lw95CL,
                 MalObsMeanLen_hi95CL=MalObsMeanLen_hi95CL,
                 MalObsMeanLen_lw95PL=MalObsMeanLen_lw95PL,
                 MalObsMeanLen_hi95PL=MalObsMeanLen_hi95PL)

  return(results)
}


#' Calculate negative log-likelihood associated with fit of growth model allowing for divergent growth between sexes or areas
#'
#' Calculates the negative log-likelihood associated with a sample of fish length-at-age data
#' and associated growth curve parameter values, for a growth model allowing for divergent
#' growth, after a certain age, between sexes or fish in different areas
#'
#' @keywords internal
#'
#' @param params c(log(c(Linf,Linf)),log(c(vbK,vbK)),tzero1,log(tdiverge))
#'
#' @return Negative log-likelihood
CalcLengthAtAge_DivergeGrowthModel <- function(params) {

  Linf = exp(params[1:2])
  if (length(params)==6) { # separate vbk
    vbK = exp(params[3:4])
    tzero1 = params[5]
    tdiverge = exp(params[6])
  }

  if (length(params)==5) { # common vbK
    vbK = rep(0,2)
    vbK[1] = exp(params[3]); vbK[2] = exp(params[3])
    tzero1 = params[4]
    tdiverge = exp(params[5])
  }

  ExpLen <- rep(NA,length(ObsAge))

  Ldiverge = Linf[1] * (1 - exp(-vbK[1] * (tdiverge - tzero1)))
  if (Ldiverge > Linf[2] - 1)  Ldiverge = Linf[2] - 1
  tzero2 = tzero1 - (1/vbK[1])*log(1-(Ldiverge/Linf[1])) +
    (1/vbK[2])*log(1-(Ldiverge/Linf[2]))

  for (j in 1:nSamples) {
    sex = ObsSex[j]
    if (sex == 1 | ObsAge[j] <= tdiverge) {
      ExpLen[j] = Linf[1] * (1 - exp(-vbK[1] * (ObsAge[j] - tzero1)))
    }
    if (sex == 2 & ObsAge[j] > tdiverge)
      ExpLen[j] = Linf[2] * (1 - exp(-vbK[2] * (ObsAge[j] - tzero2)))
  }

  return(ExpLen)
}

#' Calculate negative log-likelihood associated with fit of a growth model allowing for divergent growth
#'
#' Calculates the negative log-likelihood associated with a sample of fish length-at-age data
#' and associated growth curve parameter values, for a growth models allowing for divergent
#' growth between sexes (or areas)
#'
#' @keywords internal
#'
#' @param params c(log(c(Linf,Linf)),log(c(vbK,vbK)),tzero,log(tdiverge))
#'
#' @return Negative log-likelihood associated with growth curve fit to length-at-age data
CalcNLL_DivergeGrowthModel <- function(params) {

  nObs = ObsAge
  ExpLen=CalcLengthAtAge_DivergeGrowthModel(params)
  SqResid = ((ObsLen - ExpLen) ^ 2)
  sumSqResid = sum(SqResid)
  stdev = sqrt(sumSqResid / nObs)
  NLL = (nObs/2.) * (log(2 * pi) + 2 * log(stdev) + 1)
  results = NLL
  return(results)
}


#' Calculate negative log-likelihood associated with fit of a growth curve
#'
#' Calculates the negative log-likelihood associated with a sample of fish length-at-age data
#' and associated growth curve parameter values. Growth models include von Bertalanffy,
#' Schnute and Somers seasonal growth curve.
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero) single or combined sex, or c(log(c(Linf,Linf)),log(c(vbK,vbK)),c(tzero,tzero))
#'
#' @return Negative log-likelihood associated with growth curve fit to length-at-age data
CalcNLL_GrowthCurve <- function(params) {

  # single or separate sexes, sex recorded in data
  if (DataType == 1 | DataType == 2) {
    if (length(params) == 3 | length(params) == 6)  {
      ExpLen = CalcLengthAtAge_vonBertalanffyGrowthCurve(params) # 3, 6 or 7 params
    }
    if (length(params) == 4 | length(params) == 8) {
      ExpLen = CalcLengthAtAge_SchnuteGrowthCurve(params, t1, t2, ObsAge) # 4 or 8 params
    }
    if (length(params) == 5 | length(params) == 10) {
      ExpLen = CalcLengthAtAge_SomersSeasonalGrowthCurve(params) # 5 or 10 params
    }
  }
  # calculate expected length at age
  if (DataType == 3) { # two sexes, but sex not recorded in data
    ExpLen = CalcLengthAtAge_vonBertalanffyGrowthCurve(params) # 5, 6 or 7 params
  }

  if (DataType==1) { # lengths at age for individual fish
    if (length(params) <6) { # combined or single sex
      nObs = length(ObsAge)
      SqResid = ((ObsLen - ExpLen) ^ 2)
      sumSqResid = sum(SqResid)
      stdev = sqrt(sumSqResid / nObs)
      NLL = (nObs/2.) * (log(2 * pi) + 2 * log(stdev) + 1)
      Objfunc = NLL
    } else { # separate sexes, but sex known
      sumSqResid = 0
      nObs = 0
      for (i in 1:nSexes) {
        tempSampSize = length(which(!is.na(ObsLen[i,])))
        nObs = nObs + tempSampSize # age classes
        SqResid = ((ObsLen[i,1:tempSampSize] - ExpLen[i,1:tempSampSize]) ^ 2)
        sumSqResid = sumSqResid + sum(SqResid)
      }
      stdev = sqrt(sumSqResid / nObs)
      NLL = (nObs/2.) * (log(2 * pi) + 2 * log(stdev) + 1)
      Objfunc = NLL
    }
  }

  if (DataType==3) { # separate sexes, with sex unknown
    nObs = length(ObsLen)
    if (length(params) == 5) stdev = rep(exp(params[5]),nObs)
    if (length(params) == 6) stdev = rep(exp(params[6]),nObs)
    if (length(params) == 7) stdev = rep(exp(params[7]),nObs)
    ExpLen_1 = as.vector(ExpLen[1,])
    ExpLen_2 = as.vector(ExpLen[2,])
    LL = sum(log((dnorm(ObsLen, ExpLen_1,stdev, log=F) + dnorm(ObsLen, ExpLen_2,stdev, log=F))))
    Objfunc = -LL
  } # data type

  if (DataType==2) { # lengths at age for individual fish
    Objfunc = 0
    nObs = length(ObsAge) # age classes
    if (nSexes==1) { # combined or single sex
      for (j in 1:nObs) {
        NLL = 0.5 * log(2 * pi) + log(ObsMeanLense[j]) + ((ObsMeanLen[j] - ExpLen[j])^2) / (2 * ObsMeanLense[j]^2)
        Objfunc = Objfunc + NLL
      }
    } else { # separate sexes
      for (i in 1:nSexes) {
        for (j in 1:nObs) {
          NLL = 0.5 * log(2 * pi) + log(ObsMeanLense[i,j]) + ((ObsMeanLen[j] - ExpLen[j])^2) / (2 * ObsMeanLense[i,j]^2)
          Objfunc = Objfunc + NLL
        }
      }
    }
  }

  results = Objfunc
  return(results)
}


#' Calculate expected lengths at age from von Bertalanffy growth curve
#'
#' Calculates expected lengths at age from von Bertalanffy growth curve
#' and associated growth curve parameter values. Used to calculated expected lengths
#' for observed ages
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero) single or combined sex,
#' or c(log(c(Linf,Linf)),log(c(vbK,vbK)),c(tzero,tzero))
#' or c(log(c(Linf,Linf)),log(vbK),c(tzero,tzero))
#' or c(log(c(Linf,Linf)),log(vbK),tzero)
#'
#' @return expected lengths at age (ExpLen)
CalcLengthAtAge_vonBertalanffyGrowthCurve <- function(params) {

  # calculate expected length at age growth, for von Bertalanffy growth curve (ages of fish in sample)
  if (length(params) == 3) { # single or combined sex
    Linf = exp(params[1])
    vbK = exp(params[2])
    tzero = params[3]
    ExpLen = Linf * (1 - exp(-vbK * (ObsAge - tzero)))
  }
  if (!is.vector(ObsAge)) { # separate sexes, sexes known in data
    ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge)))
    colnames(ExpLen) = 1:length(ObsAge)
    ExpLen = as.matrix(ExpLen)
    if (length(params) == 6) {
      Linf = exp(params[1:2])
      vbK = exp(params[3:4])
      tzero = params[5:6]
      ExpLen[1,] = Linf[1] * (1 - exp(-vbK[1] * (ObsAge[1,] - tzero[1])))
      ExpLen[2,] = Linf[2] * (1 - exp(-vbK[2] * (ObsAge[2,] - tzero[2])))
    }
  }
  if (is.vector(ObsAge) & length(params)>3) {
    ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge)))
    colnames(ExpLen) = 1:length(ObsAge)
    ExpLen = as.matrix(ExpLen)
    if (DataType!=3) {
      if (length(params)== 6) {
        Linf = exp(params[1:2])
        vbK = exp(params[3:4])
        tzero = params[5:6]
        ExpLen[1,] = Linf[1] * (1 - exp(-vbK[1] * (ObsAge - tzero[1])))
        ExpLen[2,] = Linf[2] * (1 - exp(-vbK[2] * (ObsAge - tzero[2])))
      }
    }
    if (DataType==3) {
      Linf = exp(params[1:2])
      if (length(params) == 5) {
        vbK = exp(params[3])
        tzero = params[4]
        ExpLen[1,] = Linf[1] * (1 - exp(-vbK * (ObsAge - tzero)))
        ExpLen[2,] = Linf[2] * (1 - exp(-vbK * (ObsAge - tzero)))
      }
      if (length(params) == 6) {
        vbK = exp(params[3])
        tzero = params[4:5]
        ExpLen[1,] = Linf[1] * (1 - exp(-vbK * (ObsAge - tzero[1])))
        ExpLen[2,] = Linf[2] * (1 - exp(-vbK * (ObsAge - tzero[2])))
      }
      if (length(params) == 7) {
        vbK = exp(params[3:4])
        tzero = params[5:6]
        ExpLen[1,] = Linf[1] * (1 - exp(-vbK[1] * (ObsAge - tzero[1])))
        ExpLen[2,] = Linf[2] * (1 - exp(-vbK[2] * (ObsAge - tzero[2])))
      }
    }
  }

  return(ExpLen)
}

#' Calculate expected lengths at age from von Bertalanffy growth curve
#'
#' Calculates expected lengths at age from von Bertalanffy growth curve
#' and associated growth curve parameter values. Used in plotting where
#' ages are pres-specified
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero)
#' @param DataType # 1=lengths at age data for individual fish (single sex, or sexes recorded),
#' 2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
#' (two sex, sexes not recorded)
#' @param nSexes number of sexes
#' @param ObsAge observed ages
#' @param plotages specified ages for plotting
#'
#' @return specified ages (plotages) and expected lengths at ages (plotlengths)
CalcLengthAtAge_vonBertalanffyGrowthCurve2 <- function(params, DataType, nSexes, ObsAge, plotages) {
  # for plotting - calculate expected length at age growth, for von Bertalanffy growth curve (specified age range)

  if (length(params) == 3) { # single or combined sex
    if (DataType!=3) {
      Linf = exp(params[1])
      vbK = exp(params[2])
      tzero = params[3]
      plotlengths = Linf * (1 - exp(-vbK * (plotages - tzero)))
    }
  }

  if (length(params)>3) {
    if (DataType!=3) {
      if (!is.vector(ObsAge)) { # separate sexes, sexes known in data
        if (length(params) == 6) {
          plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
          colnames(plotlengths) = 1:length(plotages)
          plotlengths = as.matrix(plotlengths)
          Linf = exp(params[1:2])
          vbK = exp(params[3:4])
          tzero = params[5:6]
          plotlengths[1,] = Linf[1] * (1 - exp(-vbK[1] * (plotages - tzero[1])))
          plotlengths[2,] = Linf[2] * (1 - exp(-vbK[2] * (plotages - tzero[2])))
        }
      }

      if (is.vector(ObsAge)) {
        plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
        colnames(plotlengths) = 1:length(plotages)
        plotlengths = as.matrix(plotlengths)
        if (length(params)== 6) {
          Linf = exp(params[1:2])
          vbK = exp(params[3:4])
          tzero = params[5:6]
          plotlengths[1,] = Linf[1] * (1 - exp(-vbK[1] * (plotages - tzero[1])))
          plotlengths[2,] = Linf[2] * (1 - exp(-vbK[2] * (plotages - tzero[2])))
        }
      }
    }
  } # > 3 params

  if (DataType==3) { # sexes unknown
    plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
    colnames(plotlengths) = 1:length(plotages)
    plotlengths = as.matrix(plotlengths)
    Linf = exp(params[1:2])
    if (length(params) == 5) {
      vbK = exp(params[3])
      tzero = params[4]
      plotlengths[1,] = Linf[1] * (1 - exp(-vbK * (plotages - tzero)))
      plotlengths[2,] = Linf[2] * (1 - exp(-vbK * (plotages - tzero)))
    }
    if (length(params) == 6) {
      vbK = exp(params[3])
      tzero = params[4:5]
      plotlengths[1,] = Linf[1] * (1 - exp(-vbK * (plotages - tzero[1])))
      plotlengths[2,] = Linf[2] * (1 - exp(-vbK * (plotages - tzero[2])))
    }
    if (length(params) == 7) {
      vbK = exp(params[3:4])
      tzero = params[5:6]
      plotlengths[1,] = Linf[1] * (1 - exp(-vbK[1] * (plotages - tzero[1])))
      plotlengths[2,] = Linf[2] * (1 - exp(-vbK[2] * (plotages - tzero[2])))
    }
  }

  results = list(plotages=plotages,
                 plotlengths=plotlengths)

  return(results)

}



#' Calculate expected lengths at age from a growth model allowing for divergent growth
#'
#' Calculates expected lengths at age from von Bertalanffy growth model
#' allowing for divergent growth between sexes (or areas), used for plotting
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero, log(tdiverge))
#' @param plotages specified ages for plotting
#'
#' @return specified ages (plotages) and expected lengths at ages (plotlengths)
CalcLengthAtAge_DivergentGrowthModel2 <- function(params, plotages) {

  # for plotting - calculate expected length at age growth (for specified age range)
  plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
  colnames(plotlengths) = 1:length(plotages)
  plotlengths = as.matrix(plotlengths)
  params = as.vector(unlist(params))
  Linf = exp(params[1:2])
  if (length(params)==6) { # separate vbk
    vbK = exp(params[3:4])
    tzero1 = params[5]
    tdiverge = exp(params[6])
  }

  if (length(params)==5) { # common vbK
    vbK = rep(0,2)
    vbK[1] = exp(params[3]); vbK[2] = exp(params[3])
    tzero1 = params[4]
    tdiverge = exp(params[5])
  }

  Ldiverge = Linf[1] * (1 - exp(-vbK[1] * (tdiverge - tzero1)))
  if (Ldiverge > Linf[2] - 1)  Ldiverge = Linf[2] - 1
  tzero2 = tzero1 - (1/vbK[1])*log(1-(Ldiverge/Linf[1])) +
    (1/vbK[2])*log(1-(Ldiverge/Linf[2]))
  plotlengths[1,] = Linf[1] * (1 - exp(-vbK[1] * (plotages - tzero1)))
  plotlengths[2,] = Linf[1] * (1 - exp(-vbK[1] * (plotages - tzero1)))
  x=which(plotages > tdiverge)
  plotlengths[2,x] = Linf[2] * (1 - exp(-vbK[2] * (plotages[x] - tzero2)))

  results = list(plotages=plotages,
                 plotlengths=plotlengths)
  return(results)
}


#' Fit a von Bertalanffy growth model allowing for divergent growth between sexes or areas
#'
#' This function fits a von Bertalanffy growth curve, allowing for divergent growth between sexes
#' or areas after a certain age (tdiverge), to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb
#'
#' @keywords internal
#'
#' @param params c(log(c(Linf,Linf)),log(c(vbK,vbK)),tzero,log(tdiverge))
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsSex observed sex
#'
#' @return stored output from internal R nlminb optimisation function (nlmb)
FitDivergeGrowthModel <- function(params, ObsAge, ObsLen, ObsSex) {

  nlmb <- nlminb(params, CalcNLL_DivergeGrowthModel, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))

  results=nlmb
  return(results)
}

#' Fit a von Bertalanffy growth curve to a sample of fish length-at-age data.
#'
#' This function fits a von Bertalanffy growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero) single or combined sex, or c(log(c(Linf,Linf)),log(c(vbK,vbK)),c(tzero,tzero))
#' @param nSexes 1=single or combined sex, 2=separate sexes
#' @param DataType # 1=lengths at age data for individual fish (single sex, or sexes recorded),
#' 2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
#' (two sex, sexes not recorded)
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return stored output from internal R nlminb optimisation function (nlmb)
FitvonBertalanffyGrowthModel <- function(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  nlmb <- nlminb(params, CalcNLL_GrowthCurve, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))

  results=nlmb
  return(results)
}

#' Get statistical outputs from a fitted von Bertalanffy growth curve.
#'
#' This function fits a von Bertalanffy growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95% confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package
#'
#' @param params c(log(c(Linf,Linf)),log(c(vbK,vbK)),tzero,log(tdiverge))
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsSex observed sex (1=female, 2=male)
#'
#' @return nll, convergence, SampleSize ParamEst, PtEst_LdivergePtEst, PtEst_Maltzero, params, vcov.params, cor.params
#'
#' @examples
#' set.seed(123)
#' GrowthEqn = 4 # von Bertalanffy growth equation with divergent growth
#' nSexes = 2 # single or combined sex
#' nSamples = 500
#' MinAge = 0
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(1000,500)
#' vbK = c(0.15,0.15)
#' tzero1 = 0
#' tdiverge = 1
#' Growth_params = c(Linf, vbK, tzero1, tdiverge)
#' Ref_ages = NA
#' Growth_cv = 0.1
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' x=which(Res$ObsDat$ObsSex==1)
#' plot(Res$ObsAge[x], Res$ObsLen[x])
#' y=which(Res$ObsDat$ObsSex==2)
#' points(Res$ObsAge[y], Res$ObsLen[y],col="blue")
#' # fit growth model
#' InitLinf = c(1050,500)
#' # InitvbK = c(0.13,0.13) # separate vbK option
#' InitvbK = 0.13 # single vbK option
#' Inittzero1 = 0
#' Inittdiverge = 4
#' params = c(log(InitLinf), log(InitvbK), Inittzero1, log(Inittdiverge))
#' ObsAge = Res$ObsDat$ObsAge
#' ObsLen = Res$ObsDat$ObsLen
#' ObsSex = Res$ObsDat$ObsSex
#' FittedRes=GetDivergeGrowthModelResults(params, ObsAge, ObsLen, ObsSex)
#' FittedRes$Boot_ParamEst
#' @export
GetDivergeGrowthModelResults <- function(params, ObsAge, ObsLen, ObsSex) {

  # fit growth model
  nlmb = FitDivergeGrowthModel(params, ObsAge, ObsLen, ObsSex)

  # calculate uncertainty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(nlmb$par, CalcNLL_DivergeGrowthModel)
  params = nlmb$par
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  FemEstLinf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
  MalEstLinf = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
  if (length(params)==5) {
    EstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    FemEsttzero <- c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
    Esttdiverge <- exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), vbK=round(EstvbK,2), Femtzero=round(FemEsttzero,2),
                            MalLinf=round(MalEstLinf,1), Esttdiverge=round(Esttdiverge,2)))
  }
  if (length(params)==6) {
    FemEstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    MalEstvbK = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    FemEsttzero <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
    Esttdiverge <- exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
    ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), FemvbK=round(FemEstvbK,2), Femtzero=round(FemEsttzero,2),
                            MalLinf=round(MalEstLinf,1), MalvbK=round(MalEstvbK,2), Esttdiverge=round(Esttdiverge,2)))
  }
  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  # bootstrap confidence limits
  plotages = seq(MinAge,MaxAge,0.1)
  Res=GetConfidenceLimitsForGrowthCurve(GrowthEqn, nSexes, DataType, Ref_ages, params, vcov.params, ObsAge, ObsLen, ObsSex, plotages)
  sims.params = Res$sims.params
  Boot_ParamEst <- t(apply(sims.params, 2, quantile, probs = c(0.5, 0.025, 0.975)))
  colnames(Boot_ParamEst) = c("Boot_Est","Boot_lw_95%CL","Boot_up_95%CL")

  # tzero not log-transformed
  if (length(params)==5) {
    Boot_ParamEst[1,] = exp(Boot_ParamEst[1,]) # Linf fem
    Boot_ParamEst[2,] = exp(Boot_ParamEst[2,]) # Linf mal
    Boot_ParamEst[3,] = exp(Boot_ParamEst[3,]) # common vbK
    Boot_ParamEst[5,] = exp(Boot_ParamEst[5,]) # tdiverge
  }
  if (length(params)==6) {
    Boot_ParamEst[1,] = exp(Boot_ParamEst[1,]) # Linf fem
    Boot_ParamEst[2,] = exp(Boot_ParamEst[2,]) # Linf mal
    Boot_ParamEst[3,] = exp(Boot_ParamEst[3,]) # vbK fem
    Boot_ParamEst[4,] = exp(Boot_ParamEst[4,]) # vbK mal
    Boot_ParamEst[6,] = exp(Boot_ParamEst[6,]) # tdiverge
  }
  Boot_ParamEst = round(Boot_ParamEst,3)

  nObs = length(ObsAge)

  # store value of objective function
  nll = nlmb$objective

  # store convergence value
  convergence = nlmb$convergence

  # get point estimate for Ldiverge and tzero2
  if (length(params)==5) {
    PtEst_LdivergePtEst = ParamEst[1,1] * (1 - exp(-ParamEst[2,1] * (ParamEst[5,1] - ParamEst[2,1])))

    PtEst_Maltzero = ParamEst[3,1] - (1/ParamEst[2,1])*log(1-(PtEst_LdivergePtEst/ParamEst[1,1])) +
      (1/ParamEst[4,1])*log(1-(PtEst_LdivergePtEst/ParamEst[4,1]))
  }
  if (length(params)==6) {
    PtEst_LdivergePtEst = ParamEst[1,1] * (1 - exp(-ParamEst[2,1] * (ParamEst[6,1] - ParamEst[3,1])))

    PtEst_Maltzero = ParamEst[3,1] - (1/ParamEst[2,1])*log(1-(PtEst_LdivergePtEst/ParamEst[1,1])) +
      (1/ParamEst[4,1])*log(1-(PtEst_LdivergePtEst/ParamEst[4,1]))
  }


  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = nObs,
                 ParamEst = ParamEst,
                 Boot_ParamEst = Boot_ParamEst,
                 PtEst_LdivergePtEst = PtEst_LdivergePtEst,
                 PtEst_Maltzero = PtEst_Maltzero,
                 params = nlmb$par,
                 vcov.params = vcov.params,
                 cor.params = cor.params,
                 Fem.sim.growth.est = Res$Fem.sim.growth.est,
                 Fem.sim.growth.low = Res$Fem.sim.growth.low,
                 Fem.sim.growth.up = Res$Fem.sim.growth.up,
                 Mal.sim.growth.est = Res$Mal.sim.growth.est,
                 Mal.sim.growth.low = Res$Mal.sim.growth.low,
                 Mal.sim.growth.up = Res$Mal.sim.growth.up,
                 sim.growth.xvals = Res$plotages,
                 sims.params = Res$sims,
                 Fem.sims.curves = Res$Fem.sims.curves,
                 Mal.sims.curves = Res$Mal.sims.curves)

  return(results)

}

#' Get statistical outputs from a fitted von Bertalanffy growth curve.
#'
#' This function fits a von Bertalanffy growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95% confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package
#'
#' @param params c(log(Linf),log(vbK),tzero) single or combined sex, or c(log(c(Linf,Linf)),log(c(vbK,vbK)),c(tzero,tzero))
#' @param nSexes 1=single or combined sexes, 2=separate sexes
#' @param DataType # 1=lengths at age data for individual fish (single sex, or sexes recorded),
#' 2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
#' (two sex, sexes not recorded)
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence)
#' sample size (SampleSize), growth parameter estimates with lower and upper 95%
#' confidence limits (ParamEst), point estimates for growth parameters (params)
#' and variance-covariance matrix (vcov.params)
#'
#' @examples
#' # von Bertalanffy growth equation
#' # simulate data (ignoring mortality and selectivity effects)
#' GrowthEqn = 1 # von Bertalanffy growth equation
#' nSexes = 1 # single or combined sex
#' nSamples = 300
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' Growth_params = c(Linf, vbK, tzero)
#' Ref_ages = NA
#' Growth_cv = 0.08
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(400),log(0.3),0) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAge
#' ObsLen = Res$ObsLen
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(400),log(0.3),0) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAgeCl
#' ObsMeanLen = Res$ObsMeanLen
#' ObsMeanLense = Res$ObsMeanLense
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' nSexes = 2 # Separate sexes
#' nSamples = c(300,200) # if different, order so NAs are last for smaller sample
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(300,250)
#' vbK = c(0.3,0.3)
#' tzero = c(0,0)
#' GrowthEqn = 1
#' Growth_params = c(Linf, vbK, tzero)
#' Ref_ages = NA
#' Growth_cv = c(0.08,0.08)
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' ObsAge = as.matrix(Res$ObsAge)
#' ObsLen = as.matrix(Res$ObsLen)
#' params = c(log(c(300,300)),log(c(0.3,0.3)),c(0,0)) # log(Linf), log(k), tzero
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(c(300,300)),log(c(0.3,0.3)),c(0,0)) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAgeCl
#' FemObsMeanLen = Res$FemObsMeanLen
#' MalObsMeanLen = Res$MalObsMeanLen
#' ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' FemObsMeanLense = Res$FemObsMeanLense
#' MalObsMeanLense = Res$MalObsMeanLense
#' ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' # Simulate data for separate sexes, but fit model with 2 curves, where sex is not known
#' # 2 sexes
#' GrowthEqn = 1
#' nSamples = c(300,300)
#' nSexes = 2
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(400,500)
#' vbK = c(0.3,0.3)
#' tzero = c(0,0)
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = c(0.08,0.08)
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # plot lengths at age
#' for (i in 1:nSexes) {
#'   ObsAge=as.vector(unlist(Res$ObsAge[i,]))
#'   ObsLen=as.vector(unlist(Res$ObsLen[i,]))
#'   if (i==1) {
#'     plot(ObsAge, ObsLen, ylim=c(0,800))
#'   } else {
#'     points(ObsAge, ObsLen, col="blue")
#'   }
#' }
#' ObsAge=c(as.vector(unlist(Res$ObsAge[1,])),as.vector(unlist(Res$ObsAge[2,])))
#' ObsLen=c(as.vector(unlist(Res$ObsLen[1,])),as.vector(unlist(Res$ObsLen[2,])))
#' DataType=3  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' nSexes = 2 # separate sexes, but sex unknown
#' # note, model has 7, not 6 params, with estimated sd, to distinguish this from standard model
#' # params = c(log(c(250,800)),log(c(0.3,0.3)),c(0,0), log(20)) # log(Linf), log(k), tzero, common sd
#' # params = c(log(c(250,800)),log(0.3),c(0,0), log(20)) # log(Linf), log(k), tzero, common sd
#' params = c(log(c(250,800)),log(0.3),0, log(20)) # log(Linf), log(k), tzero, common sd
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' FittedRes$ParamEst
#' @export
GetvonBertalanffyGrowthResults <- function(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  # fit growth model
  nlmb = FitvonBertalanffyGrowthModel(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)

  # get estimates
  nlmb$objective # value of nll
  nlmb$convergence
  nlmb$par

  # calculate uncertainty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(nlmb$par, CalcNLL_GrowthCurve)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  if (nSexes==1) {
    EstLinf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstvbK = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Esttzero <- c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
    ParamEst = t(data.frame(Linf=round(EstLinf,1), vbK=round(EstvbK,2), tzero=round(Esttzero,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
    nObs = length(ObsAge)

  }
  if (nSexes==2) {
    if (!is.vector(ObsAge)) { # separate sexes, but sex known
      FemEstLinf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      MalEstLinf = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      FemEstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
      MalEstvbK = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
      FemEsttzero <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
      MalEsttzero <- c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
      ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), FemvbK=round(FemEstvbK,2), Femtzero=round(FemEsttzero,2),
                              MalLinf=round(MalEstLinf,1), MalvbK=round(MalEstvbK,2), Maltzero=round(MalEsttzero,2)))
    }
    if (is.vector(ObsAge)) { # separate sexes, with sex unknown
      FemEstLinf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
      MalEstLinf = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
      if (length(params)==5) { # common k, common tzero, sep Linf
        EstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        Esttzero <- c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
        Estsd <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
        ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), MalLinf=round(MalEstLinf,1), vbK=round(EstvbK,2),
                                tzero=round(Esttzero,2), sd=round(Estsd,2)))
      }
      if (length(params)==6) { # common k
        EstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        FemEsttzero <- c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
        MalEsttzero <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
        Estsd <- c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
        ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), MalLinf=round(MalEstLinf,1), vbK=round(EstvbK,2),
                                Femtzero=round(FemEsttzero,2), Maltzero=round(MalEsttzero,2), sd=round(Estsd,2)))
      }
      if (length(params)==7) { # separate k
        FemEstvbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        MalEstvbK = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        FemEsttzero <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
        MalEsttzero <- c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
        Estsd <- c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7])
        ParamEst = t(data.frame(FemLinf=round(FemEstLinf,1), MalLinf=round(MalEstLinf,1), FemvbK=round(FemEstvbK,2),
                                MalvbK=round(MalEstvbK,2), Femtzero=round(FemEsttzero,2), Maltzero=round(MalEsttzero,2),
                                sd=round(Estsd,2)))
      }
    }
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

    if (!is.vector(ObsAge)) {
      nObs = length(which(!is.na(ObsAge[1,]))) + length(which(!is.na(ObsAge[2,])))
    } else {
      nObs = length(ObsAge) * 2 # age classes for 2 sexes
    }
  }

  # store value of objective function
  nll = nlmb$objective

  # store convergence value
  convergence = nlmb$convergence

  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = nObs,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.params = vcov.params,
                 cor.params = cor.params)

  return(results)

}


#' Plot fish length-at-age data
#'
#' Plots either raw data (PlotCLs=FALSE) or
#' mean lengths at age with specified standard deviations for each mean length
#' (PlotCLs=TRUE)
#'
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLensd sd for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ymax maximum value for y axis
#' @param xmax maximum value for x axis
#' @param yint interval for y axis
#' @param xint interval for x axis
#' @param GraphTitle graph title
#' @param xaxis_lab label for x axis
#' @param yaxis_lab label for y axis
#' @return scatter plot of length-at-age data
#' @examples
#' # plot mean lengths at age with approx 95% CLs given sd at length
#'nSamples = 300
#'ObsAge = 1:20
#'MeanLen = 300 * (1 - exp(-0.3 * (ObsAge - 0)))
#'ObsLen = MeanLen
#'ObsMeanLensd = MeanLen*0.1
#'PlotLengthAtAgeData(ObsAge, ObsLen, ObsMeanLensd, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                    GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=TRUE)
#' # plot raw lengths at age
#'set.seed(123)
#'nSamples = 300
#'ObsAge = runif(nSamples, 1, 20)
#'MeanLen = 300 * (1 - exp(-0.3 * (ObsAge - 0)))
#'ObsLen = rnorm(nSamples, MeanLen, 0.1 * MeanLen)
#'PlotLengthAtAgeData(ObsAge, ObsLen, ObsMeanLensd=NA, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                    GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=TRUE)
#' @export
PlotLengthAtAgeData <- function(ObsAge, ObsLen, ObsMeanLensd, ymax, xmax,
                                yint, xint, GraphTitle, xaxis_lab, yaxis_lab, PlotCLs) {

  # get default axis limits and intervals
  xlims=Get_xaxis_scale(ObsAge)
  ylims=Get_yaxis_scale(ObsLen)

  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "Age, yrs"
  if (is.na(yaxis_lab)) yaxis_lab = "Length, mm"

  plot(ObsAge, ObsLen, "p", xlim=c(0,xmax), ylim=c(0,ymax), cex=0.5, pch=16,
       frame=F, xaxt = 'n', yaxt = 'n', xlab="", ylab="", main=GraphTitle)
  AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=1,  cexaxisval=NA, lwdval=NA,
                             lineval=NA, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
  mtext(xaxis_lab,las=1,side=1,line=3,cex=1.0)
  mtext(yaxis_lab,las=3,side=2,line=3,cex=1.0)
  if (PlotCLs==TRUE) {
    ObsLenlw = ObsLen - (1.96 * ObsMeanLensd)
    ObsLenup = ObsLen + (1.96 * ObsMeanLensd)
    arrows(ObsAge, ObsLenlw, ObsAge, ObsLenup,
           code=3, angle=90,length=0.02, col='black')
  }
}

#' Get confidence limits for the fitted growth curve.
#'
#' Use parametric resampling to get confidence limits for the fitted growth curve
#'
#' @param GrowthEqn 1=von Bertalanffy, 2=Schnute, 3=Somers seasonal curve, 4=divergent growth model
#' @param nSexes number of sexes
#' @param DataType 1=lengths at age data for individual fish (single sex, or sexes recorded),
#' 2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
#' (two sex, sexes not recorded)
#' @param Ref_ages Schnute growth curve reference ages (t1 and t2)
#' @param params estimated parameters
#' @param vcov.params estimated variance-covariance matrix
#' @param ObsAge observed individual fish age data (used for Growth eqn option 4)
#' @param ObsLen observed individual fish length data (used for Growth eqn option 4)
#' @param ObsSex observed individual fish sex (or area, inputted as sex) data (used for Growth eqn option 4)
#' @param plotages ages specified to produce estimates of uncertainty for estimated lengths
#'
#' @return scatter plot of length-at-age data
#'
#' @export
GetConfidenceLimitsForGrowthCurve <- function(GrowthEqn, nSexes, DataType, Ref_ages, params, vcov.params, ObsAge, ObsLen, ObsSex, plotages) {

  # store estimated parameter distributions
  set.seed(123)
  sim.growth.est=NA; sim.growth.low=NA; sim.growth.up=NA
  Fem.sim.growth.est=NA; Fem.sim.growth.low=NA; Fem.sim.growth.up=NA
  Mal.sim.growth.est=NA; Mal.sim.growth.low=NA; Mal.sim.growth.up=NA
  sims.params=NA  # required for GrowthEqn=4

  if (GrowthEqn == 1) { # von Bertalanffy
    sims = data.frame(MASS::mvrnorm(n = 1000, params, vcov.params))
    if (nSexes==1) { # single or combined sex
      GetEstLen <- function(params) {
        res=CalcLengthAtAge_vonBertalanffyGrowthCurve2(params, DataType, nSexes, ObsAge, plotages)
        return(res$plotlengths)
      }
      names(sims) = c("Linf", "vbK","tzero")
      sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLen)
    } else { # separate sexes
      GetEstLenFem <- function(params) {
        res=CalcLengthAtAge_vonBertalanffyGrowthCurve2(params, DataType, nSexes, ObsAge, plotages)
        FemEstLengths = res$plotlengths[1,]
        return(FemEstLengths)
      }
      GetEstLenMal <- function(params) {
        res=CalcLengthAtAge_vonBertalanffyGrowthCurve2(params, DataType, nSexes, ObsAge, plotages)
        MalEstLengths = res$plotlengths[2,]
        return(MalEstLengths)
      }
      if (!is.vector(ObsAge)) { # sep sexes, these known in data
        names(sims) = c("FemLinf","FemvbK","Femtzero","MalLinf","MalvbK","Maltzero")
      } else { # sep sexes, these unknown in data
        if (length(params)==5) {
          names(sims) = c("FemLinf","MalLinf","vbK","tzero","Commonsd")
        }
        if (length(params)==6) {
          names(sims) = c("FemLinf","MalLinf","vbK","Femtzero","Maltzero","Commonsd")
        }
        if (length(params)==7) { # sep sexes, these not known in data
          names(sims) = c("FemLinf","MalLinf","FemvbK","MalvbK","Femtzero","Maltzero","Commonsd")
        }
      }
      Fem.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenFem)
      Mal.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenMal)
    }
  }

  if (GrowthEqn == 2) { # Schnute
    sims = data.frame(MASS::mvrnorm(n = 1000, params, vcov.params))
    if (nSexes==1) { # single or combined sex
      GetEstLen <- function(params) {

        res=CalcLengthAtAge_SchnuteGrowthCurve2(params, nSexes, Ref_ages, plotages)
        res$plotlengths
        return(res$plotlengths)
      }
      names(sims) = c("y1", "y2","a","b")
      sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLen)
    } else { # separate sexes
      GetEstLenFem <- function(params) {
        res=CalcLengthAtAge_SchnuteGrowthCurve2(params, nSexes, Ref_ages, plotages)
        FemEstLengths = res$plotlengths[1,]
        return(FemEstLengths)
      }
      GetEstLenMal <- function(params) {
        res=CalcLengthAtAge_SchnuteGrowthCurve2(params, nSexes, Ref_ages, plotages)
        MalEstLengths = res$plotlengths[2,]
        return(MalEstLengths)
      }
      names(sims) = c("Fem_y1","Fem_y2","Fem_a","Fem_b","Mal_y1","Mal_y2","Mal_a","Mal_b")
      Fem.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenFem)
      Mal.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenMal)
    }
  }

  if (GrowthEqn == 3) { # Somers seasonal growth curve
    sims = data.frame(MASS::mvrnorm(n = 1000, params, vcov.params))
    if (nSexes==1) { # single or combined sex
      GetEstLen <- function(params) {
        res=CalcLengthAtAge_SomersSeasonalGrowthCurve2(params, nSexes, plotages)
        res$plotlengths
        return(res$plotlengths)
      }
      names(sims) = c("Linf", "vbK","tzero","tc","C")
      sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLen)
    } else { # separate sexes
      GetEstLenFem <- function(params) {
        res=CalcLengthAtAge_SomersSeasonalGrowthCurve2(params, nSexes, plotages)
        FemEstLengths = res$plotlengths[1,]
        return(FemEstLengths)
      }
      GetEstLenMal <- function(params) {
        res=CalcLengthAtAge_SomersSeasonalGrowthCurve2(params, nSexes, plotages)
        MalEstLengths = res$plotlengths[2,]
        return(MalEstLengths)
      }
      names(sims) = c("Fem_Linf","Fem_vbK","Fem_tzero","Fem_tc","Fem_C","Mal_Linf","Mal_vbK","Mal_tzero","Mal_tc","Mal_C")
      Fem.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenFem)
      Mal.sims.curves = apply(X=sims[,], MARGIN=1, FUN=GetEstLenMal)
    }
  }

  if (GrowthEqn != 4) {
    if (nSexes==1) {
      sim.growth.est = apply(sims.curves, 1, function(x) quantile(x, 0.5))
      sim.growth.low = apply(sims.curves, 1, function(x) quantile(x, 0.025))
      sim.growth.up = apply(sims.curves, 1, function(x) quantile(x, 0.975))
    }
    if (nSexes==2 | GrowthEqn==4) {
      Fem.sim.growth.est = apply(Fem.sims.curves, 1, function(x) quantile(x, 0.5))
      Fem.sim.growth.low = apply(Fem.sims.curves, 1, function(x) quantile(x, 0.025))
      Fem.sim.growth.up = apply(Fem.sims.curves, 1, function(x) quantile(x, 0.975))
      Mal.sim.growth.est = apply(Mal.sims.curves, 1, function(x) quantile(x, 0.5))
      Mal.sim.growth.low = apply(Mal.sims.curves, 1, function(x) quantile(x, 0.025))
      Mal.sim.growth.up = apply(Mal.sims.curves, 1, function(x) quantile(x, 0.975))
    }
  }

  # bootstrap resampling
  if (GrowthEqn == 4) {

    # Set number of bootstrap iterations
    OrigObsAge = ObsAge
    OrigObsLen = ObsLen
    OrigObsSex = ObsSex

    nTrials <- 200  # or higher for stable CI estimates
    sims.params <- as.data.frame(matrix(NA, nrow = nTrials, ncol = length(params)))
    if (length(params)==6) {
      names(sims.params) = c("FemLinf", "MalLinf","FemvbK","MalvbK","Femtzero","tdiverge")
    }
    if (length(params)==5) {
      names(sims.params) = c("FemLinf", "MalLinf","vbK", "Femtzero","tdiverge")
    }

    Fem.sims.curves <- matrix(NA, nrow = nTrials, ncol = length(plotages))
    Mal.sims.curves <- matrix(NA, nrow = nTrials, ncol = length(plotages))

    k=0
    for (n in 1:nTrials) {
      k=k+1
      if (k>=10) {
        k=0
        cat("Bootstrap trial", n, "of",nTrials, "\n")
        Sys.sleep(0.2)
      }

      # Resample indices
      resample_idx <- sample(1:length(ObsLen), replace = TRUE)
      ObsAge <<- OrigObsAge[resample_idx]
      ObsLen <<- OrigObsLen[resample_idx]
      ObsSex <<- OrigObsSex[resample_idx]

      # Fit the model on the bootstrap sample
      fit <- tryCatch({
        FitDivergeGrowthModel(params, ObsAge, ObsLen, ObsSex)
      }, error = function(e) {
        cat("Error on bootstrap", n, ":", e$message, "\n")
        return(NULL)
      })

      # Save parameters if fitting was successful
      if (!is.null(fit)) {
        sims.params[n, ] <- fit$par  # assuming your function returns a list with $par
      }

      params = fit$par
      res=CalcLengthAtAge_DivergentGrowthModel2(params, plotages)
      Fem.sims.curves[n,] = res$plotlengths[1,]
      Mal.sims.curves[n,] = res$plotlengths[2,]

    } # Trials

    # 95% confidence limits
    Fem.sim.growth.est = apply(Fem.sims.curves, 2, quantile, probs = 0.5, na.rm = TRUE)
    Fem.sim.growth.low = apply(Fem.sims.curves, 2, quantile, probs = 0.025, na.rm = TRUE)
    Fem.sim.growth.up = apply(Fem.sims.curves, 2, quantile, probs = 0.975, na.rm = TRUE)
    Mal.sim.growth.est = apply(Mal.sims.curves, 2, quantile, probs = 0.5, na.rm = TRUE)
    Mal.sim.growth.low = apply(Mal.sims.curves, 2, quantile, probs = 0.025, na.rm = TRUE)
    Mal.sim.growth.up = apply(Mal.sims.curves, 2, quantile, probs = 0.975, na.rm = TRUE)

  } # GrowthEqn = 4

  results = list(sims.params = sims.params,
                 sim.growth.est = sim.growth.est,
                 sim.growth.low = sim.growth.low,
                 sim.growth.up = sim.growth.up,
                 Fem.sim.growth.est = Fem.sim.growth.est,
                 Fem.sim.growth.low = Fem.sim.growth.low,
                 Fem.sim.growth.up = Fem.sim.growth.up,
                 Mal.sim.growth.est = Mal.sim.growth.est,
                 Mal.sim.growth.low = Mal.sim.growth.low,
                 Mal.sim.growth.up = Mal.sim.growth.up,
                 sim.growth.xvals = plotages)

  return(results)
}

#' Plot fitted growth curve to fish length-at-age data.
#'
#' @param DataType # 1=lengths at age data for individual fish (single sex, or sexes recorded),
#' 2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
#' (two sex, sexes not recorded)
#' @param nSexes number of sexes
#' @param GrowthEqn 1=von Bertalanffy, 2=Schnute, 3=Somers seasonal growth curve, 4=Divergent growth curve
#' @param ObsAge observed individual fish ages
#' @param ObsLen observed individual fish lengths
#' @param ObsLen observed individual fish sex (required for GrowthEqn=4), otherwise set to NA
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLensd sd for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param params estimated growth parameters
#' @param plotages specified ages for plotting
#' @param ymax maximum value for y axis
#' @param xmax maximum value for x axis
#' @param yint interval for y axis
#' @param xint interval for x axis
#' @param GraphTitle graph title
#' @param xaxis_lab label for x axis
#' @param yaxis_lab label for y axis
#' @param PlotCLs logical True=plot 95 percent confidence limits for curve
#'
#' @return fitted curve on scatter plot with length-at-age data
#' @examples
#' # von Bertalanffy growth equation
#' # simulate data (ignoring mortality and selectivity effects)
#' GrowthEqn = 1 # von Bertalanffy growth equation
#' nSexes = 1 # single or combined sex
#' nSamples = 300
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' Growth_params = c(Linf, vbK, tzero)
#' Ref_ages = NA
#' Growth_cv = 0.08
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(400),log(0.3),0) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAge
#' ObsLen = Res$ObsLen
#'
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' plotages=seq(0, MaxAge,0.1)
#' par(mfrow=c(1,1))
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen=NA,
#'                       ObsMeanLense=NA, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(400),log(0.3),0) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAgeCl
#' ObsMeanLen = Res$ObsMeanLen
#' ObsMeanLense = Res$ObsMeanLense
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' nSexes = 2 # Separate sexes
#' nSamples = c(500,500) # if different, order so NAs are last for smaller sample
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(300,250)
#' vbK = c(0.3,0.3)
#' tzero = c(0,0)
#' GrowthEqn = 1
#' Growth_params = c(Linf, vbK, tzero)
#' Ref_ages = NA
#' Growth_cv = c(0.08,0.08)
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' ObsAge = as.matrix(Res$ObsAge)
#' ObsLen = as.matrix(Res$ObsLen)
#' params = c(log(c(300,300)),log(c(0.3,0.3)),c(0,0)) # log(Linf), log(k), tzero
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' plotages=seq(0, MaxAge,0.1)
#' par(mfrow=c(1,2))
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(c(300,300)),log(c(0.3,0.3)),c(0,0)) # log(Linf), log(k), tzero
#' ObsAge = Res$ObsAgeCl
#' FemObsMeanLen = Res$FemObsMeanLen
#' MalObsMeanLen = Res$MalObsMeanLen
#' ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' FemObsMeanLense = Res$FemObsMeanLense
#' MalObsMeanLense = Res$MalObsMeanLense
#' ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Schnute growth equation - single sex length-at-age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' GrowthEqn = 2 # Schnute growth equation
#' nSexes = 1 # single or combined sex
#' nSamples = 500
#' MinAge = 0
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = c(0.5,15)
#' Growth_params = c(50,400,0.2,0.5)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' ObsAge=Res$ObsAge
#' ObsLen=Res$ObsLen
#' par(mfrow=c(1,1))
#' plot(ObsAge,ObsLen)
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = log(c(80,400,0.2,0.5))
#' t1=0.5
#' t2=15
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' FittedRes$ParamEst
#' FittedRes$cor.params
#' plotages=seq(0, MaxAge,0.1)
#' par(mfrow=c(1,1))
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen=NA,
#'                       ObsMeanLense=NA, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Schnute growth equation - single sex mean length at age data
#' # simulate data (ignoring mortality and selectivity effects)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' nSamples = 2000
#' MinAge = 0
#' MaxAge = 3
#' AgeStep = 0.1
#' Ref_ages = c(0.1,2.5)
#' Growth_params = c(20,200,1.0,0.5)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' ObsAge=Res$ObsAge
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLense=Res$ObsMeanLense
#' plot(Res$ObsAgeCl,ObsMeanLen)
#' params = log(c(25,180,0.9,0.5))
#' t1=0.1
#' t2=2.5
#' ObsAge=Res$ObsAgeCl
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLense=Res$ObsMeanLense
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' FittedRes$ParamEst
#' FittedRes$cor.params
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Schnute growth equation - 2 sexes, length-at-age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' GrowthEqn = 2 # Schnute equation
#' nSexes = 2 # separate sexes
#' nSamples = c(500,500)
#' MinAge = 1
#' MaxAge = 25
#' AgeStep = 1
#' Ref_ages = c(1,15)
#' Growth_params = c(100,120,500,550,0.2,0.2,0.3,0.4)
#' Growth_cv = c(0.1,0.1)
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # Fit Schnute growth curve to length at age dat
#' ObsAge=as.matrix(Res$ObsAge)
#' ObsLen=as.matrix(Res$ObsLen)
#' par(mfrow=c(1,1))
#' plot(ObsAge[1,],ObsLen[1,])
#' points(ObsAge[2,],ObsLen[2,],col="blue")
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = log(c(90,110,420,430,0.3,0.3,0.6,0.6))
#' t1=1
#' t2=15
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' FittedRes$ParamEst
#' FittedRes$cor.params
#' par(mfrow=c(1,2))
#' plotages=seq(0, MaxAge,1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen=NA,
#'                       ObsMeanLense=NA, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Fit Schnute growth curve - 2 sexes mean length at age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' nSamples = c(2000,2000)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' MinAge = 0.1
#' MaxAge = 3
#' AgeStep = 0.1
#' ObsAge = runif(nSamples, MinAge, MaxAge)
#' Growth_params = c(20,20,180,200,1.0,1.0,0.5,0.5)
#' Ref_ages = c(0.1,2.5)
#' Growth_cv = c(0.1,0.1)
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' par(mfrow=c(1,1))
#' plot(Res$ObsAgeCl, Res$FemObsMeanLen)
#' points(Res$ObsAgeCl, Res$MalObsMeanLen,col="blue")
#' t1=0.1
#' t2=2.5
#' params = log(c(20,30,190,210,0.9,1.1,0.5,0.5))
#' ObsAge=Res$ObsAgeCl
#' FemObsMeanLen = Res$FemObsMeanLen
#' MalObsMeanLen = Res$MalObsMeanLen
#' ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' FemObsMeanLense = Res$FemObsMeanLense
#' MalObsMeanLense = Res$MalObsMeanLense
#' ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' FittedRes$ParamEst
#' FittedRes$cor.params
#' par(mfrow=c(1,2))
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Seasonal growth curve
#' # simulate data (ignoring mortality and selectivity effects)
#' GrowthEqn=3 # Seasonal growth curve
#' nSexes=1 # single or combined sex
#' nSamples = 1000
#' MinAge = 0
#' MaxAge = 3
#' AgeStep = 1/12
#' Linf = 150
#' vbK = 1
#' tzero = 0
#' tc = 0.25
#' C = 0.8
#' Ref_ages=NA
#' Growth_params=c(Linf,vbK,tzero,tc,C)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' DataType = 1
#' nSexes=1
#' ObsAge = Res$ObsAge
#' ObsLen = Res$ObsLen
#' params = c(log(160),log(0.8),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' # fit seasonal growth curve to mean length at age data
#' FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' par(mfrow=c(1,1))
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen=NA, ObsMeanLense=NA,
#'                       params, Ref_ages, plotages, ymax=200, xmax=NA, yint=50, xint=NA, GraphTitle=NA,
#'                       xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(160),log(0.8),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' ObsAge=Res$ObsAgeCl
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLense=Res$ObsMeanLense
#' FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' GrowthEqn=3 # Seasonal growth curve
#' nSexes=2 # separate sexes
#' nSamples = c(500,500)
#' MinAge = 0
#' MaxAge = 3
#' AgeStep = 1/12
#' Linf = c(150,160)
#' vbK = c(1,1)
#' tzero = c(0,0)
#' tc = c(0.25,0.25)
#' C = c(0.8,0.8)
#' Ref_ages=NA
#' Growth_params=c(Linf,vbK,tzero,tc,C)
#' Growth_cv = c(0.1,0.1)
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' DataType = 1
#' nSexes=2
#' ObsAge=as.matrix(Res$ObsAge)
#' ObsLen=as.matrix(Res$ObsLen)
#' params = c(log(c(160,150)),log(c(0.8,0.8)),c(0,0),c(0.25,0.25),c(0.8,0.8)) #Linf, vbK, tzero, tc, C
#' # fit seasonal growth curve to mean length at age data
#' FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' par(mfrow=c(1,2))
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen=NA,
#'                       ObsMeanLense=NA, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = c(log(c(160,150)),log(c(0.8,0.8)),c(0,0),c(0.25,0.25),c(0.8,0.8)) #Linf, vbK, tzero, tc, C
#' ObsAge=Res$ObsAgeCl
#' FemObsMeanLen = Res$FemObsMeanLen
#' MalObsMeanLen = Res$MalObsMeanLen
#' ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' FemObsMeanLense = Res$FemObsMeanLense
#' MalObsMeanLense = Res$MalObsMeanLense
#' ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' plotages=seq(0, MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # simulate data allowing for divergent growth
#' set.seed(123)
#' GrowthEqn = 4 # von Bertalanffy growth equation with divergent growth
#' nSexes = 2 # single or combined sex
#' nSamples = 500
#' MinAge = 0
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(1000,500)
#' vbK = c(0.15,0.15)
#' tzero1 = 0
#' tdiverge = 1
#' Growth_params = c(Linf, vbK, tzero1, tdiverge)
#' Ref_ages = NA
#' Growth_cv = 0.1
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' x=which(Res$ObsDat$ObsSex==1)
#' plot(Res$ObsAge[x], Res$ObsLen[x])
#' y=which(Res$ObsDat$ObsSex==2)
#' points(Res$ObsAge[y], Res$ObsLen[y],col="blue")
#' # fit growth model
#' InitLinf = c(1050,500)
#' # InitvbK = c(0.13,0.13) # separate vbK option
#' InitvbK = 0.13 # single vbK option
#' Inittzero1 = 0
#' Inittdiverge = 4
#' params = c(log(InitLinf), log(InitvbK), Inittzero1, log(Inittdiverge))
#' ObsAge = Res$ObsDat$ObsAge
#' ObsLen = Res$ObsDat$ObsLen
#' ObsSex = Res$ObsDat$ObsSex
#' FittedRes=GetDivergeGrowthModelResults(params, ObsAge, ObsLen, ObsSex)
#' FittedRes$Boot_ParamEst
#' DataType=1
#' plotages = seq(MinAge,MaxAge,0.1)
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' # Simulate data for separate sexes, fit separate von Bertalanffy growth curves, where sex is not specified in data
#' # 2 sexes
#' GrowthEqn = 1
#' nSamples = c(300,300)
#' nSexes = 2
#' MinAge = 1
#' MaxAge = 20
#' AgeStep = 1
#' Linf = c(400,600)
#' vbK = c(0.3,0.3)
#' tzero = c(0,0)
#' Growth_params = c(Linf,vbK,tzero)
#' Growth_cv = c(0.08,0.08)
#' Res = SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # plot lengths at age
#' for (i in 1:nSexes) {
#'   ObsAge=as.vector(unlist(Res$ObsAge[i,]))
#'   ObsLen=as.vector(unlist(Res$ObsLen[i,]))
#'   if (i==1) {
#'     plot(ObsAge, ObsLen, ylim=c(0,800))
#'   } else {
#'     points(ObsAge, ObsLen, col="blue")
#'   }
#' }
#' ObsAge=c(as.vector(unlist(Res$ObsAge[1,])),as.vector(unlist(Res$ObsAge[2,])))
#' ObsLen=c(as.vector(unlist(Res$ObsLen[1,])),as.vector(unlist(Res$ObsLen[2,])))
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' nSexes = 2 # separate sexes, but sex unknown
#' # note, model has 7, not 6 params, with estimated sd, to distinguish this from standard model
#' params = c(log(c(250,800)),log(c(0.3,0.3)),c(0,0), log(20)) # log(Linf), log(k), tzero, common sd
#' FittedRes=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' plotages=seq(0, MaxAge,0.1)
#' par(mfrow=c(1,1),mar=c(4,3,2,2))
#' PlotFittedGrowthCurve(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex=NA, ObsMeanLen,
#'                       ObsMeanLense, params, Ref_ages, plotages, ymax=NA, xmax=NA, yint=NA, xint=NA,
#'                       GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T, FittedRes)
#' @export
PlotFittedGrowthCurve <- function(DataType, nSexes, GrowthEqn, ObsAge, ObsLen, ObsSex, ObsMeanLen,
                                  ObsMeanLense, params, Ref_ages, plotages, ymax, xmax, yint, xint,
                                  GraphTitle, xaxis_lab, yaxis_lab, PlotCLs, FittedRes) {


  # get default axis limits and intervals
  if(is.vector(ObsAge)) {
    xlims=Get_xaxis_scale(ObsAge)
    if (DataType==1 | DataType==3) ylims=Get_yaxis_scale(ObsLen)
    if (DataType==2) ylims=Get_yaxis_scale(ObsMeanLen)
  }
  if(is.matrix(ObsAge)) {
    newObsAge = as.vector(c(ObsAge[1,],ObsAge[2,]))
    newObsAge = newObsAge[which(!is.na(newObsAge))]
    xlims=Get_xaxis_scale(newObsAge)
  }
  if (DataType==1) {
    if(is.matrix(ObsLen)) {
      newObsLen = as.vector(c(ObsLen[1,],ObsLen[2,]))
      newObsLen = newObsLen[which(!is.na(newObsLen))]
      ylims=Get_yaxis_scale(newObsLen)
    }
  }
  if (DataType==2) {
    if(is.matrix(ObsMeanLen)) {
      newObsLen = as.vector(c(ObsMeanLen[1,],ObsMeanLen[2,]))
      newObsLen = newObsLen[which(!is.na(newObsLen))]
      ylims=Get_yaxis_scale(newObsLen)
    }
  }

  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "Age, yrs"
  if (is.na(yaxis_lab)) yaxis_lab = "Length, mm"

  # fit growth curve and get results, including variance-covariance matrix
  if (GrowthEqn == 1) { # von Bertalanffy
    if (is.list(FittedRes)) {     # if model already fitted, can input results rather than refit
      Res =  FittedRes
    } else {
      Res=GetvonBertalanffyGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
    }
  }
  if (GrowthEqn == 2) { # Schnute
    if (is.list(FittedRes)) {
      Res =  FittedRes
    } else {
      Res=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
    }
  }
  if (GrowthEqn == 3) { # seasonal growth
    if (is.list(FittedRes)) {
      Res =  FittedRes
    } else {
      Res=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
    }
  }
  if (GrowthEqn == 4) { # divergent growth
    if (is.list(FittedRes)) {
      Res =  FittedRes
    } else {
      Res=GetDivergeGrowthModelResults(params, ObsAge, ObsLen, ObsSex)
    }
  }

  params = Res$params
  vcov.params = Res$vcov.params
  if (GrowthEqn != 4) { # not divergent growth (as already have results from Get function)
    Res=GetConfidenceLimitsForGrowthCurve(GrowthEqn, nSexes, DataType, Ref_ages, params, vcov.params,
                                          ObsAge, ObsLen, ObsSex, plotages)
  }

  if (GrowthEqn == 1) res = CalcLengthAtAge_vonBertalanffyGrowthCurve2(params, DataType, nSexes, ObsAge, plotages) # von Bertalanffy
  if (GrowthEqn == 2) res = CalcLengthAtAge_SchnuteGrowthCurve2(params, nSexes, Ref_ages, plotages) # Schnute
  if (GrowthEqn == 3) res = CalcLengthAtAge_SomersSeasonalGrowthCurve2(params, nSexes, plotages) # Seasonal growth
  if (GrowthEqn == 4) res = CalcLengthAtAge_DivergentGrowthModel2(params, plotages) # Divergent growth

  if (DataType == 1 | DataType==3) FishLen = ObsLen
  if (DataType == 2) FishLen = ObsMeanLen

  if (nSexes==1 | GrowthEqn == 4 | DataType==3) { # single or combined sex
    plot(ObsAge, FishLen, "p", xlim=c(0,xmax), ylim=c(0,ymax), cex=0.4, pch=16,
         frame=F, xaxt = 'n', yaxt = 'n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=1.2,  cexaxisval=NA, lwdval=NA,
                               lineval=0, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2)
    mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2)

    if (GrowthEqn==3) {
      # single sex, seasonal growth curve
      lines(res$plotages, res$plotlengths, lwd=2)
    }
    if (GrowthEqn!=3) {
      if (GrowthEqn == 4 | length(params)==5 | length(params)==6 | length(params)==7) {
        # sex specific but sex not recorded or divergent growth model
        if (GrowthEqn !=3) {
          lines(res$plotages, res$plotlengths[1,], lwd=2)
          lines(res$plotages, res$plotlengths[2,], lwd=2)
        }
      } else {
        # single sex
        lines(res$plotages, res$plotlengths, lwd=2)
      }
    }

    if (PlotCLs == TRUE) { # plot confidence limits
      # single sex, seasonal growth curve
      if (GrowthEqn==3) {
        x = c(res$plotages,rev(res$plotages)) # using shading for 95% CLs
        y = c(Res$sim.growth.low, rev(Res$sim.growth.up))
        polygon(x,y,col="light grey",border=NA)
        lines(res$plotages, Res$sim.growth.est, "l", lty="solid")
        lines(res$plotages, Res$sim.growth.low, "l", lty="dotted")
        lines(res$plotages, Res$sim.growth.up, "l", lty="dotted")
        points(ObsAge, FishLen, cex=0.4, pch=16)
      }
      # sex specific but sex not recorded or divergent growth model
      if (GrowthEqn!=3) {
        if (GrowthEqn == 4 | length(params)==5 | length(params)==6 | length(params)==7) {
          x = c(res$plotages,rev(res$plotages)) # using shading for 95% CLs
          y = c(Res$Fem.sim.growth.low, rev(Res$Fem.sim.growth.up))
          polygon(x,y,col="pink",border=NA)
          y = c(Res$Mal.sim.growth.low, rev(Res$Mal.sim.growth.up))
          polygon(x,y,col="light blue",border=NA)
          lines(res$plotages, Res$Fem.sim.growth.est, "l", lty="solid", col="red")
          lines(res$plotages, Res$Fem.sim.growth.low, "l", lty="dotted", col="red")
          lines(res$plotages, Res$Fem.sim.growth.up, "l", lty="dotted", col="red")
          lines(res$plotages, Res$Mal.sim.growth.est, "l", lty="solid", col="blue")
          lines(res$plotages, Res$Mal.sim.growth.low, "l", lty="dotted", col="blue")
          lines(res$plotages, Res$Mal.sim.growth.up, "l", lty="dotted", col="blue")
          points(ObsAge, FishLen, cex=0.4, pch=16)

        } else {
          # single sex
          x = c(res$plotages,rev(res$plotages)) # using shading for 95% CLs
          y = c(Res$sim.growth.low, rev(Res$sim.growth.up))
          polygon(x,y,col="light grey",border=NA)
          lines(res$plotages, Res$sim.growth.est, "l", lty="solid")
          lines(res$plotages, Res$sim.growth.low, "l", lty="dotted")
          lines(res$plotages, Res$sim.growth.up, "l", lty="dotted")
          points(ObsAge, FishLen, cex=0.4, pch=16)
        }
      }
    }
    if (DataType == 2) {
      ObsLenlw = ObsMeanLen - (1.96 * ObsMeanLense)
      ObsLenup = ObsMeanLen + (1.96 * ObsMeanLense)
      arrows(ObsAge, ObsLenlw, ObsAge, ObsLenup, code=3, angle=90,length=0.02, col='black')
    }
  } else { # separate sexes
    # females
    if (DataType == 1) AgesToPlot = ObsAge[1,]
    if (DataType == 2 | DataType==3) AgesToPlot = ObsAge
    plot(AgesToPlot, FishLen[1,], "p", xlim=c(0,xmax), ylim=c(0,ymax), cex=0.4, pch=16,
         frame=F, xaxt = 'n', yaxt = 'n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=1.2,  cexaxisval=NA, lwdval=NA,
                               lineval=0, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2)
    mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2)
    lines(res$plotages, res$plotlengths[1,], lwd=2)
    legend("topleft", pch=-1, legend="Females", bty='n', cex=0.8,lwd=-1)
    if (PlotCLs == TRUE) { # plot confidence limits
      x = c(res$plotages,rev(res$plotages)) # using shading for 95% CLs
      y = c(Res$Fem.sim.growth.low, rev(Res$Fem.sim.growth.up))
      polygon(x,y,col="pink",border=NA)
      lines(res$plotages, Res$Fem.sim.growth.est, "l", lty="solid", col="red")
      lines(res$plotages, Res$Fem.sim.growth.low, "l", lty="dotted", col="red")
      lines(res$plotages, Res$Fem.sim.growth.up, "l", lty="dotted", col="red")
      points(AgesToPlot, FishLen[1,], cex=0.4, pch=16)
    }
    if (DataType == 2) {
      ObsLenlw = ObsMeanLen[1,] - (1.96 * ObsMeanLense[1,])
      ObsLenup = ObsMeanLen[1,] + (1.96 * ObsMeanLense[1,])
      arrows(ObsAge, ObsLenlw, ObsAge, ObsLenup,
             code=3, angle=90,length=0.02, col='black')
    }
    # males
    if (DataType == 1) AgesToPlot = ObsAge[2,]
    if (DataType == 2 | DataType==3) AgesToPlot = ObsAge
    plot(AgesToPlot, FishLen[2,], "p", xlim=c(0,xmax), ylim=c(0,ymax), cex=0.4, pch=16,
         frame=F, xaxt = 'n', yaxt = 'n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax, yint, cexval=1.2,  cexaxisval=NA, lwdval=NA,
                               lineval=0, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2)
    mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2)
    lines(res$plotages, res$plotlengths[2,], lwd=2)
    legend("topleft", pch=-1, legend="Males", bty='n', cex=0.8,lwd=-1)
    if (PlotCLs == TRUE) {  # plot confidence limits
      x = c(res$plotages,rev(res$plotages)) # using shading for 95% CLs
      y = c(Res$Mal.sim.growth.low, rev(Res$Mal.sim.growth.up))
      polygon(x,y,col="light blue",border=NA)
      lines(res$plotages, Res$Mal.sim.growth.est, "l", lty="solid", col="blue")
      lines(res$plotages, Res$Mal.sim.growth.low, "l", lty="dotted", col="blue")
      lines(res$plotages, Res$Mal.sim.growth.up, "l", lty="dotted", col="blue")
      points(AgesToPlot, FishLen[2,], cex=0.4, pch=16)
    }
    if (DataType == 2) {
      ObsLenlw = ObsMeanLen[2,] - (1.96 * ObsMeanLense[2,])
      ObsLenup = ObsMeanLen[2,] + (1.96 * ObsMeanLense[2,])
      arrows(ObsAge, ObsLenlw, ObsAge, ObsLenup,
             code=3, angle=90,length=0.02, col='black')
    }
  }
}


#' Calculate estimated length at age from Schnute growth curve
#'
#' Calculates estimated length at age from Schnute growth curve, given the age, growth parameters and
#' reference ages
#'
#' Function requires two reference ages, t1 and t2 (stored in memory in R).
#' @param Age specified age
#' @param t1 first reference age
#' @param t2 second reference age
#' @param y1 length at t1
#' @param y2 length at t2
#' @param a growth curve parameter
#' @param a growth curve parameter
#'
#' @return Estimated length at a specified age
#'
#' @examples
#' t1=0.5; t2=10
#' y1=50; y2=250; a=0; b=3
#' Age=2
#' SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
#' @export
SchnuteGrowthfunction <- function (Age, t1, t2, y1, y2, a, b) {
  # Schnute's versatile growth model

  # Reference:
  #    Schnute, J. T. and Richards, L. J.  (1990).  A unified aproach to the
  #    analysis of fish growth, maturity and survivorship data.  Can. J.
  #    Fish. Aquat. Sci. 47: 24-40.

  # If the age lies below the theoretical age at which the length is zero,
  # the estimated length is assumed to be zero.

  # robustify routine, prevent y2 being less than y1
  if (y2 < y1+2) {
    y2 = y1+2
    # cat("SchnuteGrowthfunction: y1", y1,"y2",y2,'\n')
  }

  # ' Determine which equation is to be used
  if (a == 0) {
    if (b == 0) {
      # Eqn(18)
      y = y1 * exp(log(y2 / y1) * (Age - t1) / (t2 - t1))

      if (y < 0) {
        y = 0
      } # y < 0
    } else {          #(b == 0)
      # Eqn(17)
      # First, let's work out tzero
      #cat("Age", Age, "y1",y1,"y2",y2,"b",b,'\n')
      tzero = t1 - (y1 ^ b) * (t2 - t1) / (y2 ^ b - y1 ^ b)

      if (Age < tzero) {
        y = 0
      } else {
        v = (y1 ^ b + (y2 ^ b - y1 ^ b) * (Age - t1) / (t2 - t1))
        y = v ^ (1 / b)
      }
    }
  }# a == 0

  if (a != 0) {
    if (b == 0) {
      # Eqn(16)
      y = y1 * exp(log(y2 / y1) * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
      if (y < 0) {
        y = 0 }
    } else {
      # Eqn(15)
      # First. let's work out tzero
      if (1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b) <= 0) {
        tzero = t1 - log(1E-4) / a
      } else {
        tzero = t1 - log(1 + (y1 ^ b) * (1 - exp(-a * (t2 - t1))) / (y2 ^ b - y1 ^ b)) / a
      }
      if (is.nan(tzero)) {
        cat("SchnuteGrowthfunction: Problem calculating tzero",'\n')
      }
      # cat("Age",Age, "tzero",tzero, '\n')
      if (Age < tzero) {
        y = 0
      } else {
        v = (y1 ^ b + (y2 ^ b - y1 ^ b)
             * (1 - exp(-a * (Age - t1))) / (1 - exp(-a * (t2 - t1))))
        y = v ^ (1 / b)
      } # else
    } # b == 0
  }  # a != 0

  return(y)
} # end function

#' Calculate expected lengths at age from Schnute growth curve
#'
#' Calculates expected lengths at age from Schnute growth curve
#' and associated growth curve parameter values. Used to calculated expected lengths
#' for observed ages
#'
#' @keywords internal
#'
#' @param params c(log(y1),log(y2),a,b) Schnute parameters
#' @param t1 reference age 1
#' @param t2 reference age 2
#' @param ObsAge specified ages
#'
#' @return expected lengths at age (ExpLen)
#' @export
CalcLengthAtAge_SchnuteGrowthCurve <- function(params, t1, t2, ObsAge) {

  if (length(params) == 4) nSexes = 1
  if (length(params) == 8) nSexes = 2

  # calculate expected length at age growth, for Schnute growth curve (ages of fish in sample)
  if (nSexes==1) { # single or combined sex
    y1 = exp(params[1])
    y2 = exp(params[2])
    a = exp(params[3])
    b = exp(params[4])
    nObs <- length(ObsAge)
    ExpLen = rep(NA,nObs)
    for (i in 1:nObs) {
      Age <- ObsAge[i]
      ExpLen[i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
    }
  } else { # separate sexes
    y1 = exp(params[1:2])
    y2 = exp(params[3:4])
    a = exp(params[5:6])
    b = exp(params[7:8])

    if (DataType==1) { # lengths at age
      ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge[1,])))
      colnames(ExpLen) = 1:length(ObsAge[1,])
      ExpLen = as.matrix(ExpLen)
      for (i in 1:nSexes) {
        nObs <- length(which(!is.na(ObsAge[i,])))
        for (j in 1:nObs) {
          Age <- ObsAge[i,j]
          ExpLen[i,j] = SchnuteGrowthfunction(Age, t1, t2, y1[i], y2[i], a[i], b[i])
        }
      }
    } else { # mean lengths at ages
      ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge)))
      colnames(ExpLen) = 1:length(ObsAge)
      ExpLen = as.matrix(ExpLen)
      for (i in 1:nSexes) {
        nObs <- length(ObsAge)
        for (j in 1:nObs) {
          Age <- ObsAge[j]
          ExpLen[i,j] = SchnuteGrowthfunction(Age, t1, t2, y1[i], y2[i], a[i], b[i])
        }
      }
    }
  }

  return(ExpLen)
}


#' Calculate expected lengths at age from Schnute growth curve
#'
#' Calculates expected lengths at age from Schnute growth curve
#' and associated growth curve parameter values. Used to calculated expected lengths
#' at specified ages for plotting.
#'
#' @keywords internal
#'
#' @param params c(log(y1),log(y2),a,b) Schnute parameters
#' @param nSexes number of sexes
#' @param Ref_ages Schnute growth curve reference ages (t1 and t2)
#' @param ObsAge specified ages
#' @param plotages specified ages for plotting
#'
#' @return specified ages (plotages) and expected lengths at ages (plotlengths)
CalcLengthAtAge_SchnuteGrowthCurve2 <- function(params, nSexes, Ref_ages, plotages) {

  nObs = length(plotages)
  t1 = Ref_ages[1]
  t2 = Ref_ages[2]

  if (length(params) == 4) nSexes = 1
  if (length(params) == 8) nSexes = 2

  if (nSexes==1) { # single or combined sex
    y1 = exp(params[1])
    y2 = exp(params[2])
    a = exp(params[3])
    b = exp(params[4])
    plotlengths = rep(NA,nObs)
    for (i in 1:nObs) {
      Age <- plotages[i]
      plotlengths[i] = SchnuteGrowthfunction(Age, t1, t2, y1, y2, a, b)
    }
  } else { # separate sexes
    y1 = exp(params[1:2])
    y2 = exp(params[3:4])
    a = exp(params[5:6])
    b = exp(params[7:8])
    plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
    colnames(plotlengths) = 1:length(plotages)
    plotlengths = as.matrix(plotlengths)
    for (i in 1:nObs) {
      Age <- plotages[i]
      plotlengths[1,i] = SchnuteGrowthfunction(Age, t1, t2, y1[1], y2[1], a[1], b[1])
      plotlengths[2,i] = SchnuteGrowthfunction(Age, t1, t2, y1[2], y2[2], a[2], b[2])
    }
  }

  results = list(plotages=plotages,
                plotlengths=plotlengths)

  return(results)
}


#' Fit a Schnute growth curve to a sample of fish length-at-age data.
#'
#' This function fits a Schnute growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' @keywords internal
#'
#' @param params c(log(y1),log(y2),a,b) Schnute parameters
#' @param DataType 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' @param nSexes number of sexes
#' @param t1 first reference age
#' @param t2 second reference age
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitSchnuteGrowthModel <- function(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  nlmb <- nlminb(params, CalcNLL_GrowthCurve, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))

  results=nlmb
  return(results)
}

#' Get outputs from a fitted Schnute growth curve.
#'
#' This function fits a Schnute growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95 percent confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package.
#'
#' @param params c(log(y1),log(y2),a,b) Schnute parameters
#' @param DataType 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' @param nSexes number of sexes
#' @param t1 first reference age
#' @param t2 second reference age
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence),
#' sample size (SampleSize), growth parameter estimates with lower and upper 95 percent confidence limits
#' (ParamEst), point estimates for growth parameters (params), variance-covariance matrix (vcov.params)
#'
#' @examples
#' # Schnute growth equation - single sex length-at-age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' GrowthEqn = 2 # Schnute growth equation
#' nSexes = 1 # single or combined sex
#' nSamples = 500
#' MinAge = 0
#' MaxAge = 20
#' AgeStep = 1
#' Ref_ages = c(0.5,15)
#' Growth_params = c(50,400,0.2,0.5)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' ObsAge=Res$ObsAge
#' ObsLen=Res$ObsLen
#' par(mfrow=c(1,1))
#' plot(ObsAge,ObsLen)
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = log(c(80,400,0.2,0.5))
#' t1=0.5
#' t2=15
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' FittedRes$ParamEst
#' # Schnute growth equation - single sex mean length at age data
#' # simulate data (ignoring mortality and selectivity effects)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' nSamples = 2000
#' MinAge = 0
#' MaxAge = 3
#' AgeStep = 0.1
#' Ref_ages = c(0.1,2.5)
#' Growth_params = c(20,200,1.0,0.5)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model and plot
#' ObsAge=Res$ObsAge
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLense=Res$ObsMeanLense
#' plot(Res$ObsAgeCl,ObsMeanLen)
#' params = log(c(25,180,0.9,0.5))
#' t1=0.1
#' t2=2.5
#' ObsAge=Res$ObsAgeCl
#' ObsMeanLen=Res$ObsMeanLen
#' ObsMeanLense=Res$ObsMeanLense
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' FittedRes$ParamEst
#' # Schnute growth equation - 2 sexes, length-at-age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' GrowthEqn = 2 # Schnute equation
#' nSexes = 2 # separate sexes
#' nSamples = c(500,500)
#' MinAge = 1
#' MaxAge = 25
#' AgeStep = 1
#' Ref_ages = c(1,15)
#' Growth_params = c(100,120,500,550,0.2,0.2,0.3,0.4)
#' Growth_cv = c(0.1,0.1)
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # Fit Schnute growth curve to length at age dat
#' ObsAge=as.matrix(Res$ObsAge)
#' ObsLen=as.matrix(Res$ObsLen)
#' par(mfrow=c(1,1))
#' plot(ObsAge[1,],ObsLen[1,])
#' points(ObsAge[2,],ObsLen[2,],col="blue")
#' DataType=1  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' params = log(c(90,110,420,430,0.3,0.3,0.6,0.6))
#' t1=1
#' t2=15
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen=NA, ObsMeanLense=NA)
#' FittedRes$ParamEst
#' # Fit Schnute growth curve - 2 sexes mean length at age data
#' # simulate data (ignoring mortality and selectivity effects)
#' set.seed(123)
#' nSamples = c(2000,2000)
#' DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' MinAge = 0.1
#' MaxAge = 3
#' AgeStep = 0.1
#' ObsAge = runif(nSamples, MinAge, MaxAge)
#' Growth_params = c(20,20,180,200,1.0,1.0,0.5,0.5)
#' Ref_ages = c(0.1,2.5)
#' Growth_cv = c(0.1,0.1)
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' par(mfrow=c(1,1))
#' plot(Res$ObsAgeCl, Res$FemObsMeanLen)
#' points(Res$ObsAgeCl, Res$MalObsMeanLen,col="blue")
#' t1=0.1
#' t2=2.5
#' params = log(c(20,30,190,210,0.9,1.1,0.5,0.5))
#' ObsAge=Res$ObsAgeCl
#' FemObsMeanLen = Res$FemObsMeanLen
#' MalObsMeanLen = Res$MalObsMeanLen
#' ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' FemObsMeanLense = Res$FemObsMeanLense
#' MalObsMeanLense = Res$MalObsMeanLense
#' ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' FittedRes=GetSchnuteGrowthResults(params, nSexes, DataType, t1, t2, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' FittedRes$ParamEst
#' @export
GetSchnuteGrowthResults <- function(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  # fit growth model
  GrowthEqn = 2
  nlmb = FitSchnuteGrowthModel(params, nSexes, DataType, t1, t2, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)

  # get estimates
  nlmb$objective # value of nll
  nlmb$convergence
  nlmb$par

  # calculate uncertainty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(nlmb$par, CalcNLL_GrowthCurve)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlations
  cor.params = temp %*% vcov.params %*% temp

  if (nSexes==1) {
    y1 = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    y2 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    a <- exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    b <- exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    ParamEst = t(data.frame(y1=round(y1,1), y2=round(y2,1), a=round(a,3), b=round(b,3)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
    SampleSize = length(ObsAge)

  } else {
    Femy1 = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    Maly1 = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Femy2 = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    Maly2 = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    Fema <- exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    Mala <- exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
    Femb <- exp(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7]))
    Malb <- exp(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8]))
    ParamEst = t(data.frame(Fem_y1=round(Femy1,1), Fem_y2=round(Femy2,2), Fem_a=round(Fema,2), Fem_b=round(Femb,2),
                            Mal_y1=round(Maly1,1), Mal_y2=round(Maly2,2), Mal_a=round(Mala,2), Mal_b=round(Malb,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

    if (DataType==1) {
      Fem_nObs = length(which(!is.na(ObsAge[1,])))
      Mal_nObs = length(which(!is.na(ObsAge[2,])))
      SampleSize = data.frame(Fem_nObs=Fem_nObs, Mal_nObs=Mal_nObs)
    } else {
      SampleSize = length(ObsAge) * 2 # age classes for 2 sexes
    }
  }

  # store value of objective function
  nll = nlmb$objective

  # store convergence value
  convergence = nlmb$convergence

  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.params = vcov.params,
                 cor.params = cor.params)

  return(results)

}

#' Get expected lengths at age from a seasonal growth curve
#'
#' Returns expected lengths at age from Somers (1988) seasonal growth curve, given parameter values
#' and ages
#'
#' @param params c(log(Linf), log(vbK), tzero, tc, C)
#'
#' @return Returns expected lengths at age from Somers (1988) seasonal growth curve
#'
#' @examples
#' # Generate length at age data for individual fish
#' DataType = 1
#' set.seed(123)
#' nSamples = 300
#' ObsAge = runif(nSamples, 0.1, 3)
#' params = c(log(150),log(1),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' MeanLen = CalcLengthAtAge_SomersSeasonalGrowthCurve(params)
#' ObsLen = rnorm(nSamples, MeanLen, 0.1*MeanLen)
#' plot(ObsAge, ObsLen)
#' params = c(log(150),log(1),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' res=CalcLengthAtAge_SomersSeasonalGrowthCurve(params)
#' points(ObsAge, res)
#' # Generate monthly mean length (with error) data, from mixture analysis
#' DataType = 2
#' set.seed(123)
#' ObsAge = seq(0.2, 3, 1/12)
#' nAges=length(ObsAge)
#' params = c(log(150),log(1),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' ObsMeanLen = CalcLengthAtAge_SomersSeasonalGrowthCurve(params)
#' ObsMeanLensd = rnorm(nAges,5,1)
#' ObsLenlw = ObsMeanLen - (1.96 * ObsMeanLensd)
#' ObsLenup = ObsMeanLen + (1.96 * ObsMeanLensd)
#' ymax=1.1*max(ObsLenup)
#' plot(ObsAge, ObsMeanLen, ylim=c(0,ymax))
#' arrows(ObsAge, ObsLenlw, ObsAge, ObsLenup,
#'        code=3, angle=90,length=0.02, col='black')
#' params = c(log(150),log(1),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' res=CalcLengthAtAge_SomersSeasonalGrowthCurve(params)
#' lines(ObsAge,res)
#' @export
CalcLengthAtAge_SomersSeasonalGrowthCurve <- function(params) {

  if (length(params) == 5) nSexes = 1
  if (length(params) == 10) nSexes = 2

  if (nSexes==1) {
    Linf = exp(params[1])
    vbK = exp(params[2])
    tzero = params[3]
    tc = params[4]
    C = params[5]
    S_age_t = sin(2 * pi * (ObsAge - tc))
    S_tzero = sin(2 * pi * (tzero - tc))
    ExpLen = Linf * (1 - exp(-vbK * (ObsAge - tzero +
                                       (C / (2 * pi)) * (S_age_t - S_tzero))))
  } else { # separate sexes
    Linf = exp(params[1:2])
    vbK = exp(params[3:4])
    tzero = params[5:6]
    tc = params[7:8]
    C = params[9:10]
    if (DataType==1) { # lengths at age
      ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge[1,])))
      colnames(ExpLen) = 1:length(ObsAge[1,])
      ExpLen = as.matrix(ExpLen)
      for (i in 1:nSexes) {
        nObs = length(which(!is.na(ObsAge[i,]))) # allow for different sample sizes
        AgesForAnalysis = as.vector(unlist(ObsAge[i,1:nObs]))
        S_tzero = sin(2 * pi * (tzero[i] - tc[i]))
        S_age_t = sin(2 * pi * (AgesForAnalysis - tc[i]))
        ExpLen[i,1:nObs] = Linf[i] *(1 - exp(-vbK[i] * (AgesForAnalysis - tzero[i] +
                                              (C[i] / (2 * pi)) * (S_age_t - S_tzero))))
      }
    } else { # mean lengths at ages
      ExpLen = data.frame(matrix(nrow=2, ncol=length(ObsAge)))
      colnames(ExpLen) = 1:length(ObsAge)
      ExpLen = as.matrix(ExpLen)
      for (i in 1:nSexes) {
        S_tzero = sin(2 * pi * (tzero[i] - tc[i]))
        S_age_t = sin(2 * pi * (ObsAge - tc[i]))
        ExpLen[i,] = Linf[i] *(1 - exp(-vbK[i] * (ObsAge - tzero[i] +
                                                    (C[i] / (2 * pi)) * (S_age_t - S_tzero))))
      }
    }
  }

  return(ExpLen)
}


#' Fit a seasonal growth curve to a sample of fish length-at-age data
#'
#' This function fits a Somers (1988) growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb
#'
#' @keywords internal
#'
#' @param params Linf, vbK, tzero, tc, C
#' @param nSexes number of sexes
#' @param DataType 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitSeasonalGrowthModel <- function(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  nlmb <- nlminb(params, CalcNLL_GrowthCurve, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))

  results=nlmb
  return(results)
}

#' Get outputs from a fitted seasonal growth curve.
#'
#' This function fits Somers (1988) seasonal growth curve to a sample of fish length-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95% confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package.
#'
#' @param params c(log(Linf), log(vbK), tzero, tc, C)
#' @param nSexes number of sexes
#' @param DataType 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' @param ObsAge observed ages
#' @param ObsLen observed lengths
#' @param ObsMeanLen mean lengths for specified ages (i.e. as estimated from mixture analysis)
#' @param ObsMeanLense se for mean estimated mean lengths for specified ages (i.e. as estimated from mixture analysis)
#'
#' @return negative log-likelihood (nll)
#' nlminb convergence diagnostic (convergence)
#' sample size (SampleSize)
#' growth parameter estimates with lower and upper 95% confidence limits (ParamEst)
#' point estimates for growth parameters (params)
#' variance-covariance matrix (vcov.params)
#' @examples
#' # Seasonal growth curve
#' # simulate data (ignoring mortality and selectivity effects)
#' GrowthEqn=3 # Seasonal growth curve
#' nSexes=1 # single or combined sex
#' nSamples = 1000
#' MinAge = 0
#' MaxAge = 3
#' AgeStep = 1/12
#' Linf = 150
#' vbK = 1
#' tzero = 0
#' tc = 0.25
#' C = 0.8
#' Ref_ages=NA
#' Growth_params=c(Linf,vbK,tzero,tc,C)
#' Growth_cv = 0.1
#' Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # fit growth model
#' DataType = 1
#' nSexes=1
#' ObsAge = Res$ObsAge
#' ObsLen = Res$ObsLen
#' ObsMeanLen=NA
#' ObsMeanLense=NA
#' params = c(log(160),log(0.8),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' # fit seasonal growth curve to mean length at age data
#' FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
#' # DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' # params = c(log(160),log(0.8),0,0.25,0.8) #Linf, vbK, tzero, tc, C
#' # ObsAge=Res$ObsAgeCl
#' # ObsMeanLen=Res$ObsMeanLen
#' # ObsMeanLense=Res$ObsMeanLense
#' # FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
#' # GrowthEqn=3 # Seasonal growth curve
#' # nSexes=2 # separate sexes
#' # nSamples = c(500,500)
#' # MinAge = 0
#' # MaxAge = 3
#' # AgeStep = 1/12
#' # Linf = c(150,160)
#' # vbK = c(1,1)
#' # tzero = c(0,0)
#' # tc = c(0.25,0.25)
#' # C = c(0.8,0.8)
#' # Ref_ages=NA
#' # Growth_params=c(Linf,vbK,tzero,tc,C)
#' # Growth_cv = c(0.1,0.1)
#' # Res=SimulateLengthAtAgeData(GrowthEqn, nSamples, nSexes, MinAge, MaxAge, AgeStep, Ref_ages, Growth_params, Growth_cv)
#' # # fit growth model
#' # DataType = 1
#' # nSexes=2
#' # ObsAge=as.matrix(Res$ObsAge)
#' # ObsLen=as.matrix(Res$ObsLen)
#' # ObsMeanLen=NA
#' # ObsMeanLense=NA
#' # params = c(log(c(160,150)),log(c(0.8,0.8)),c(0,0),c(0.25,0.25),c(0.8,0.8)) #Linf, vbK, tzero, tc, C
#' # # fit seasonal growth curve to mean length at age data
#' # FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)
#' # DataType=2  # 1=lengths at age data for individual fish, 2=mean length at age and sd data from mixture analysis
#' # params = c(log(c(160,150)),log(c(0.8,0.8)),c(0,0),c(0.25,0.25),c(0.8,0.8)) #Linf, vbK, tzero, tc, C
#' # ObsAge=Res$ObsAgeCl
#' # FemObsMeanLen = Res$FemObsMeanLen
#' # MalObsMeanLen = Res$MalObsMeanLen
#' # ObsMeanLen = as.matrix(t(data.frame(FemObsMeanLen=FemObsMeanLen,MalObsMeanLen=MalObsMeanLen)))
#' # FemObsMeanLense = Res$FemObsMeanLense
#' # MalObsMeanLense = Res$MalObsMeanLense
#' # ObsMeanLense = as.matrix(t(data.frame(FemObsMeanLense=FemObsMeanLense,MalObsMeanLense=MalObsMeanLense)))
#' # FittedRes=GetSeasonalGrowthResults(params, nSexes, DataType, ObsAge, ObsLen=NA, ObsMeanLen, ObsMeanLense)
#' @export
GetSeasonalGrowthResults <- function(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense) {

  # fit growth model
  GrowthEqn = 3
  nlmb = FitSeasonalGrowthModel(params, nSexes, DataType, ObsAge, ObsLen, ObsMeanLen, ObsMeanLense)

  # get estimates
  nlmb$objective # value of nll
  nlmb$convergence
  nlmb$par

  # calculate uncertianty for parameter estimates by getting variance-covariance matrix,
  # from fitted model, to get standard errors
  hess.out = optimHess(nlmb$par, CalcNLL_GrowthCurve)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  if (nSexes==1) {
    EstLinf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstvbK = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Esttzero <- c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
    Esttc <- c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
    EstC <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])

    # store results in data frame
    ParamEst = t(data.frame(Linf=round(EstLinf,1), vbK=round(EstvbK,2), tzero=round(Esttzero,2), tc=round(Esttc,2), C=round(EstC,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
    SampleSize = length(ObsAge)

  } else {

    Fem_Linf = c(exp(nlmb$par[1]), exp(nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    Mal_Linf = c(exp(nlmb$par[2]), exp(nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Fem_vbK = c(exp(nlmb$par[3]), exp(nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    Mal_vbK = c(exp(nlmb$par[4]), exp(nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    Fem_tzero <- c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
    Mal_tzero <- c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
    Fem_tc <- c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7])
    Mal_tc <- c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8])
    Fem_C <- c(nlmb$par[9], nlmb$par[9] + c(-1.96, 1.96) * ses[9])
    Mal_C <- c(nlmb$par[10], nlmb$par[10] + c(-1.96, 1.96) * ses[10])
    ParamEst = t(data.frame(Fem_Linf=round(Fem_Linf,1), Fem_vbK=round(Fem_vbK,2), Fem_tzero=round(Fem_tzero,2), Fem_tc=round(Fem_tc,2), Fem_C=round(Fem_C,2),
                            Mal_Linf=round(Mal_Linf,1), Mal_vbK=round(Mal_vbK,2), Mal_tzero=round(Mal_tzero,2), Mal_tc=round(Mal_tc,2), Mal_C=round(Mal_C,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

    if (DataType==1) {
      Fem_nObs = length(which(!is.na(ObsAge[1,])))
      Mal_nObs = length(which(!is.na(ObsAge[2,])))
      SampleSize = data.frame(Fem_nObs=Fem_nObs, Mal_nObs=Mal_nObs)
    } else {
      SampleSize = length(ObsAge) * 2 # age classes for 2 sexes
    }
  }


  # store sample size
  SampleSize = length(ObsAge)

  # store value of objective function
  nll = nlmb$objective

  # store convergence value
  convergence = nlmb$convergence

  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 SampleSize = SampleSize,
                 ParamEst = ParamEst,
                 params = nlmb$par,
                 vcov.params = vcov.params,
                 cor.params = cor.params)

  return(results)

}

#' Calculate expected lengths at age from seasonal growth curve
#'
#' Calculates expected lengths at age from Somers seasonal growth curve
#' and associated growth curve parameter values. Used in plotting where
#' ages are pres-specified
#'
#' @keywords internal
#'
#' @param params c(log(Linf),log(vbK),tzero,tc,C)
#' @param params number of sexes
#' @param plotages specified ages for plotting
#'
#' @return specified ages (plotages) and expected lengths at age (plotlengths)
CalcLengthAtAge_SomersSeasonalGrowthCurve2 <- function(params, nSexes, plotages) {
  # for plotting - calculate expected length at age growth, for seasonal growth curve (specified age range)

  # define set of ages over range of data for plotting
  nObs = length(plotages)

  if (length(params) == 5) nSexes = 1
  if (length(params) == 10) nSexes = 2

  if (nSexes == 1) {
    Linf = exp(params[1])
    vbK = exp(params[2])
    tzero = params[3]
    tc = params[4]
    C = params[5]
    S_age_t = rep(NA, nObs)
    plotlengths = rep(NA, nObs)

    S_tzero = sin(2 * pi * (tzero - tc))
    for (j in 1:nObs) {
      S_age_t[j] = sin(2 * pi * (plotages[j] - tc))
      plotlengths[j] = Linf *(1 - exp(-vbK * (plotages[j] - tzero +
                                                (C / (2 * pi)) * (S_age_t[j] - S_tzero))))
    }

  } else { # separate sexes
    Linf = exp(params[1:2])
    vbK = exp(params[3:4])
    tzero = params[5:6]
    tc = params[7:8]
    C = params[9:10]
    plotlengths = data.frame(matrix(nrow=2, ncol=length(plotages)))
    colnames(plotlengths) = 1:length(plotages)
    plotlengths = as.matrix(plotlengths)
    S_age_t = plotlengths
    for (i in 1:nSexes) {
      S_tzero = sin(2 * pi * (tzero[i] - tc[i]))
      for (j in 1:nObs) {
        S_age_t[i,j] = sin(2 * pi * (plotages[j] - tc[i]))
        plotlengths[i,j] = Linf[i] *(1 - exp(-vbK[i] * (plotages[j] - tzero[i] +
                                                  (C[i] / (2 * pi)) * (S_age_t[i,j] - S_tzero))))
      }
    }
  }

  results = list(plotages=plotages,
                plotlengths=plotlengths)

  return(results)
}


#****************************
# GROWTH - tag-recapture data
#****************************

#' Simulate some tag-recapture data
#'
#' Simulate tag-recapture data with varying times at liberty,
#' to which a growth curve may be fitted
#'
#' @param GrowthCrvChoice 1=double logistic model, 2=Gaussian function, 3=von Bertalanffy, 4=Gompertz
#' @param nstep number of numerical integration steps (higher number increases accuracy but reduces program speed)
#' @param nobs number of observations
#' @param MaxLen maximum length
#' @param params log(c(L50_1, L95_1, L50_2, L95_2)) double logistic model, or
#' log(c(Gaussian_A, Gaussian_u, Gaussian_sd)) Gaussian function, or log(c(vb_Linf, vb_K)) von Bertalanffy, or
#' log(c(Gomp_Linf, Gomp_G)) Gompertz
#'
#' @return simulated tag-recapture data for individuals animals, including random a set of initial lengths (Obs_Initlen),
#' final lengths (Obs_Finlen) and durations at liberty (Obs_delta_t), and the 'true' relationships (for two
#' durations at liberty), using for generating the random tag-recapture data.
#'
#' @examples
#' # Gausian
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' Gaussian_A = 0.94
#' Gaussian_u = 55.33
#' Gaussian_sd = 53.17
#' StandDev = 10
#' GrowthCrvChoice = 2 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(Gaussian_A, Gaussian_u, Gaussian_sd, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' #
#' # von Bertalanffy
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 140
#' vb_K = 0.2
#' Stdev = 3
#' GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, vb_K, Stdev))
#' nobs = 1000
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' #
#' # Gompertz
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 172
#' g = 0.35
#' Stdev = 5
#' GrowthCrvChoice = 4 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, g, Stdev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' #
#' # double logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50_1 = 46.48
#' L95_1 = 121.31
#' L50_2 = 121.91
#' L95_2 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 1 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(L50_1, L95_1, L50_2, L95_2, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' #
#' # inverse logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50 = 121.91
#' L95 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 5 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve,
#' # 4 = Gompertz growth curve, 5=inverse logistic
#' params = log(c(L50, L95, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' @export
SimulateTagRecaptureData <- function(GrowthCrvChoice, nstep, nobs, MaxLen, params) {

  if (GrowthCrvChoice == 1) StandDev = exp(params[6]) # double logistic
  if (GrowthCrvChoice == 2) StandDev = exp(params[4]) # Gaussian function
  if (GrowthCrvChoice == 3) StandDev = exp(params[3]) # von Bertalanffy
  if (GrowthCrvChoice == 4) StandDev = exp(params[3]) # Gompertz

  EstLenAtRelAge = rep(0, nobs)
  Exp_Finlen = rep(NA, nobs)

  # generate random initial lengths
  rand_Initlen = round(runif(nobs, 0.05*MaxLen, 0.95* MaxLen),0)
  # specify times at liberty
  rand_delta_t = c(rep(365, nobs/2),rep(2*365, nobs/2))

  StartAge = 0
  Obs_delta_t = rand_delta_t
  Obs_Initlen = rand_Initlen
  CalculationStage = 1
  LenPrevIntAge=NA
  for (j in 1:nobs) {
    Exp_Finlen[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge,
                             Obs_delta_t, Obs_Initlen)
  }

  Obs_Finlen = rnorm(length(Exp_Finlen),Exp_Finlen, StandDev)

  delta_t_line = rep(365,MaxLen)
  Initlen_line = 1:MaxLen
  Exp_Finlen2 = rep(NA,MaxLen)
  LenPrevIntAge=NA
  for (j in 1:MaxLen) {
    Exp_Finlen2[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge,
                              delta_t_line, Initlen_line)
  }

  delta_t_line2 = rep(2*365,MaxLen)
  Initlen_line2 = 1:MaxLen
  Exp_Finlen3 = rep(NA,MaxLen)
  LenPrevIntAge=NA
  for (j in 1:MaxLen) {
    Exp_Finlen3[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge,
                              delta_t_line2, Initlen_line2)
  }

  results = list(Obs_delta_t = Obs_delta_t,
                 Obs_Initlen = Obs_Initlen,
                 Obs_Finlen = Obs_Finlen,
                 delta_t_line = delta_t_line,
                 Initlen_line = Initlen_line,
                 Exp_Finlen2 = Exp_Finlen2,
                 delta_t_line2 = delta_t_line2,
                 Initlen_line2 = Initlen_line2,
                 Exp_Finlen3 = Exp_Finlen3)

  return(results)

}

#' Calculate negative log-likelihood for growth model fitted to tag-recapture data.
#'
#' Calculates the negative log-likelihood associated with a sample of animal tag-recapture data
#' and associated parameters of the growth model. Four alternative growth models are currently
#' available, i.e. von Bertalanffy curve, Gompertz curve, Gaussian function curve, Double logistic
#' function curve. Each describe the expected length increment between initial capture, measurement
#' and marking, and final capture and measurement.
#' Function requires data for individual animals on size at initial capture, size at final capture,
#' and time at liberty. Obs_Initlen, Obs_Finlen, Obs_delta_t (Numeric Vectors)
#' (stored in memory in R).
#'
#' @keywords internal
#'
#' This function (with parameter inputs) can be passed into R optimisation routines (e.g. nlminb).
#' @param  Model 1: L50_1, L95_1, L50_2, L95_2, Max_increment (double logistic model)
# ' Model 2: Gaussian_A, Gaussian_u, Gaussian_sd (Gaussian function)
#'  Model 3: vb_Linf, vb_K (von Bertalanffy)
#'  Model 4: Gomp_Linf, Gomp_G
#'  Model 5: L50, L95, Max_increment
#'  Model 6: MinGrowth, MaxGrowth, RateOfGrowthChange
#' @return Negative-log likelihood associated with growth curve fit to tag-recapture data
CalcNLL_TaggingGrowthModel <- function(params) {

  if (GrowthCrvChoice == 1) StandDev = exp(params[6]) # double logistic
  if (GrowthCrvChoice == 2) StandDev = exp(params[4]) # Gaussian
  if (GrowthCrvChoice == 3) StandDev = exp(params[3]) # von Bert
  if (GrowthCrvChoice == 4) StandDev = exp(params[3]) # Gompertz
  if (GrowthCrvChoice == 5) StandDev = exp(params[4]) # Inverse logistic
  if (GrowthCrvChoice == 6) StandDev = exp(params[3]) # Exponential decay

  EstDeltaL = TaggingGrowthModelNLLCalcs_Rcpp(params, nobs, GrowthCrvChoice, nstep, Obs_delta_t, Obs_Initlen);

  NLL = -sum((dnorm(Obs_Finlen - Obs_Initlen, EstDeltaL, StandDev, log = TRUE)))

  return(NLL)
}

#' Fit a tagging growth model to tag-recapture data
#'
#' This function fits a tagging growth model to a sample of fish tag-recapture data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' @keywords internal
#' Obs_delta_t, Obs_Initlen, Obs_Finlen (Numeric Vectors)
#'
#' @param params log(c(L50_1, L95_1, L50_2, L95_2, Max_increment)) double logistic model, or
#' log(c(Gaussian_A, Gaussian_u, Gaussian_sd)) Gaussian function, or log(c(vb_Linf, vb_K)) von Bertalanffy, or
#' log(c(Gomp_Linf, Gomp_G)) Gompertz
#' @param nstep number of numerical integration steps (higher number increases accuracy but reduces program speed)
#' @param Obs_delta_t observed durations at liberty for individual animals
#' @param Obs_Initlen observed initial lengths for individual animals
#' @param Obs_Finlen observed final lengths for individual animals
#' @param nobs number of observations
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitTaggingGrowthModel <- function(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, nobs) {

  nlmb <- nlminb(params, CalcNLL_TaggingGrowthModel, gradient = NULL,
                 hessian = TRUE,  control=list(trace=1, eval.max=1000, iter.max=1000))

  results=nlmb
  return(results)
}

#' Get statistical outputs from a fitted tagging growth model.
#'
#' This function fits a model to fish tag-recapture data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95 percent confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package.
#'
#' @param params log(c(L50_1, L95_1, L50_2, L95_2, Max_increment, sd)) double logistic model, or
#' log(c(Gaussian_A, Gaussian_u, Gaussian_sd, sd)) Gaussian function, or log(c(vb_Linf, vb_K, sd)),
#' von Bertalanffy, log(c(Gomp_Linf, Gomp_G, sd)) Gompertz or log(c(L50, L95, Max_increment, sd)) Inverse logistic
#' @param nstep number of numerical integration steps (higher number increases accuracy but reduces program speed)
#' @param Obs_delta_t observed durations at liberty for individual animals
#' @param Obs_Initlen observed initial lengths for individual animals
#' @param Obs_Finlen observed final lengths for individual animals
#' @param MaxAge specified maximum age for species
#' @param MaxLen maximum length to consider for analysis
#' @param nobs number of observations
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence),
#' growth parameter estimates with lower and upper 95 percent confidence limits (ParamEst),
#' point estimates for growth parameters (params), variance-covariance matrix (vcov.params),
#' observed differences between final and initial lengths (Obs_delta_Length), observed final
#' lengths (Est_Final_Length), specified set of initial lengths for plotting model results
#' (plot_Init_Length), estimated annual growth given initial length, for plotting model results
#' (plot_Length_Inc), integer ages for plotting model results (plot_IntAge), estimated lengths
#' at integer ages for plotting (plot_EstLenAtIntAge)
#'
#' @examples
#' # Gausian
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' Gaussian_A = 0.94
#' Gaussian_u = 55.33
#' Gaussian_sd = 53.17
#' StandDev = 10
#' GrowthCrvChoice = 2 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(Gaussian_A, Gaussian_u, Gaussian_sd, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(0.1, 60, 50, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' #
#' # von Bertalanffy
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 140
#' vb_K = 0.2
#' Stdev = 3
#' GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, vb_K, Stdev))
#' nobs = 1000
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' # fit model to simulated data
#' params = log(c(160, 0.3, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' MaxAge = 40
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' #
#' # Gompertz
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 172
#' g = 0.35
#' Stdev = 5
#' GrowthCrvChoice = 4 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, g, Stdev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' # fit model to simulated data
#' params = log(c(160, 0.3, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' MaxAge = 20
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' #
#' # double logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50_1 = 46.48
#' L95_1 = 121.31
#' L50_2 = 121.91
#' L95_2 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 1 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(L50_1, L95_1, L50_2, L95_2, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(60, 110, 120, 130, 0.1, 5))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' #
#' # inverse logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50 = 121.91
#' L95 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 5 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve,
#' # 4 = Gompertz growth curve, 5=inverse logistic
#' params = log(c(L50, L95, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(120, 130, 0.1, 5))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' @export
#'
GetTaggingGrowthModelResults <- function(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs) {

  # fit maturity curve
  nlmb = FitTaggingGrowthModel(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, nobs)

  nll = nlmb$objective # store value of objective function
  convergence = nlmb$convergence  # store convergence value

  # get variance-covariance matrix, from fitted model, to get standard errors
  params = nlmb$par
  hess.out = optimHess(nlmb$par, CalcNLL_TaggingGrowthModel)
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # get asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp


  # Calculate the terms that are used to produce the estimate of the derivative
  # of length
  if (GrowthCrvChoice == 1)   { # double logistic

    # L50_1, L95_1, L50_2, L95_2, Max_increment (double logistic model)

    # calculate 95 percent confidence limits
    EstL50_1 = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstL95_1 =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstL50_2 =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstL95_2 =  exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
    Max_increment =  exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
    EstStandDev =  exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))

    # store results in data frame
    ParamEst = t(data.frame(L50_1=round(EstL50_1,2), L95_1=round(EstL95_1,2),
                            L50_2=round(EstL50_2,2), L95_2=round(EstL95_2,2),
                            Max_increment=round(Max_increment,2), StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  }

  if (GrowthCrvChoice == 2)   { # Gaussian function

    # Gaussian_A, Gaussian_u, Gaussian_sd (Gaussian function)

    EstGaussian_A =  exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstGaussian_u =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstGaussian_sd =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstStandDev =  exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))

    # store results in data frame
    ParamEst = t(data.frame(Gaussian_A=round(EstGaussian_A,2), Gaussian_u=round(EstGaussian_u,2),
                            Gaussian_sd=round(EstGaussian_sd,2),StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  }

  if (GrowthCrvChoice == 3)   { # von Bertalanffy growth curve. delta_t is
    # converted t to decimal years (i.e., 1/52)

    # vb_Linf, vb_K (von Bertalanffy)
    Estvb_Linf =  exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    Estvb_K =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstStandDev =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))

    # store results in data frame
    ParamEst = t(data.frame(vb_Linf=round(Estvb_Linf,2), vb_K=round(Estvb_K,2),
                            StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  }

  if (GrowthCrvChoice == 4)   { # Gompertz growth curve. delta_t is
    # converted to decimal years (i.e., 1/52)

    EstGomp_Linf =  exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstGomp_G =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstStandDev =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))

    # store results in data frame
    ParamEst = t(data.frame(Gomp_Linf=round(EstGomp_Linf,2), Gomp_G=round(EstGomp_G,2),
                            StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  }

  if (GrowthCrvChoice == 5)   { # Inverse logistic growth curve.

    EstL50 =  exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    EstL95 =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    Max_increment =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
    EstStandDev =  exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))

    # store results in data frame
    ParamEst = t(data.frame(EstL50=round(EstL50,2), EstL95=round(EstL95,2),
                            Max_increment=round(Max_increment,2), StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
  }

  if (GrowthCrvChoice == 6)   { # Exponential decay model

    MaxGrowth =  exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
    RateOfGrowthChange =  exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
    EstStandDev =  exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))

    # store results in data frame
    ParamEst = t(data.frame(EstMaxGrowth=round(MaxGrowth,2),
                            RateOfGrowthChange=round(RateOfGrowthChange,2),
                            StandDev=round(EstStandDev,2)))
    colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")
  }

  # deltaL vs initial length
  CalculationStage = 1
  LenPrevIntAge = NA
  EstLenAtRelAge = rep(0, nobs)
  for (j in 1:nobs) {
    EstLenAtRelAge[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge=0,
                                 Obs_delta_t, Obs_Initlen)
  }

  # estimate of annual growth increment vs specified initial length
  CalculationStage = 2
  LenInc = rep(0,MaxLen)
  for (j in 1:MaxLen)  {
    LenInc[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge=0,
                         Obs_delta_t, Obs_Initlen) - j
  }

  # estimate of length at integer ages
  CalculationStage = 3
  EstLenAtIntAge <- rep(0, MaxAge)
  for (j in seq(1,MaxAge+1,1)) {
    if (j==1) {
      EstLenAtIntAge[j] = 0 # age zero
      if (GrowthCrvChoice == 4) {
        EstLenAtIntAge[j] = 0.1 # age zero
      }
    } else {
      LenPrevIntAge = EstLenAtIntAge[j-1]
      EstLenAtIntAge[j] = LenAtAge_Rcpp(j, params, GrowthCrvChoice, nstep, CalculationStage, LenPrevIntAge, StartAge=0,
                                   Obs_delta_t, Obs_Initlen)
    }
  }

  Obs_delta_Length = Obs_Finlen-Obs_Initlen
  Est_delta_Length = EstLenAtRelAge-Obs_Initlen
  Est_Final_Length = EstLenAtRelAge
  plot_Init_Length = 1:MaxLen
  plot_Length_Inc = LenInc
  plot_IntAge = 0:MaxAge
  plot_EstLenAtIntAge = EstLenAtIntAge

  # store all results as a list object
  results = list(nll = nll,
                 convergence = convergence,
                 ParamEst = ParamEst,
                 params = params,
                 vcov.params = vcov.params,
                 cor.params = cor.params,
                 Obs_delta_Length=Obs_delta_Length,
                 Est_delta_Length=Est_delta_Length,
                 Est_Final_Length=Est_Final_Length,
                 plot_Init_Length=plot_Init_Length,
                 plot_Length_Inc=plot_Length_Inc,
                 plot_IntAge=plot_IntAge,
                 plot_EstLenAtIntAge=plot_EstLenAtIntAge)
  return(results)

}


#' Plot fitted tagging growth model and data.
#'
#' This function fits a model to fish tag-recapture data
#' by minimising the negative log-likelihood associated with the parameters and data,
#' using nlminb. It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95% confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package. Results are then plotted.
#'
#' params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, FittedRes
#'
#' @param params log(c(L50_1, L95_1, L50_2, L95_2)) double logistic model, or
#' log(c(Gaussian_A, Gaussian_u, Gaussian_sd)) Gaussian function, or log(c(vb_Linf, vb_K)) von Bertalanffy, or
#' log(c(Gomp_Linf, Gomp_G)) Gompertz
#' @param nstep number of numerical integration steps (higher number increases accuracy but reduces program speed)
#' @param Obs_delta_t observed durations at liberty for individual animals
#' @param Obs_Initlen observed initial lengths for individual animals
#' @param Obs_Finlen observed final lengths for individual animals
#' @param MaxAge specified maximum age for species
#' @param MaxLen maximum length to consider for analysis
#' @param nobs number of observation
#' @param FittedRes outputs of fitted model from nlminb (can be set to NA or saved values from the
#' GetTaggingGrowthModelResults function, to speed up plotting
#' @param PlotOpt 0=all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' 5=annual len inc vs init len, 6=len vs age
#'
#' @return various plots showing model fit to tag-recapture data, also including residual plots, and simulated length at age curve
#'
#' @examples
#' # Gausian
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' Gaussian_A = 0.94
#' Gaussian_u = 55.33
#' Gaussian_sd = 53.17
#' StandDev = 10
#' GrowthCrvChoice = 2 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(Gaussian_A, Gaussian_u, Gaussian_sd, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(0.1, 60, 50, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' params = FittedRes$params
#' PlotOpt = 0 # all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' # 5=annual len inc vs init len, 6=len vs age
#' PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt)
#' #
#' # von Bertalanffy
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 140
#' vb_K = 0.2
#' Stdev = 3
#' GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, vb_K, Stdev))
#' nobs = 1000
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' # fit model to simulated data
#' params = log(c(160, 0.3, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' MaxAge = 40
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' params = FittedRes$params
#' PlotOpt = 0 # all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' # 5=annual len inc vs init len, 6=len vs age
#' PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt)
#' #
#' # Gompertz
#' # Simulate data
#' set.seed(123)
#' nstep = 100 # number of steps for numerical integration
#' MaxLen = 300
#' vb_Linf = 172
#' g = 0.35
#' Stdev = 5
#' GrowthCrvChoice = 4 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(vb_Linf, g, Stdev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Initlen_line, res$Exp_Finlen2, col="blue", "l")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' lines(res$Obs_Initlen, res$Obs_Finlen, "p")
#' # fit model to simulated data
#' params = log(c(160, 0.3, 8))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' MaxAge = 20
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' params = FittedRes$params
#' PlotOpt = 0 # all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' # 5=annual len inc vs init len, 6=len vs age
#' PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt)
#' #
#' # double logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50_1 = 46.48
#' L95_1 = 121.31
#' L50_2 = 121.91
#' L95_2 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 1 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
#' params = log(c(L50_1, L95_1, L50_2, L95_2, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(60, 110, 120, 130, 0.1, 5))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' params = FittedRes$params
#' PlotOpt = 0 # all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' # 5=annual len inc vs init len, 6=len vs age
#' PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt)
#' #
#' # inverse logistic
#' # Simulate data
#' set.seed(123)
#' nstep = 50 # number of steps for numerical integration
#' MaxLen = 300
#' L50 = 121.91
#' L95 = 171.84
#' a = 0.113
#' StandDev = 5
#' GrowthCrvChoice = 5 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve,
#' # 4 = Gompertz growth curve, 5=inverse logistic
#' params = log(c(L50, L95, a, StandDev))
#' nobs = 200
#' res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
#' plot(res$Obs_Initlen, res$Obs_Finlen, pch=16, cex=0.6)
#' lines(res$Initlen_line, res$Exp_Finlen2, col="blue")
#' lines(res$Initlen_line2, res$Exp_Finlen3, col="blue")
#' MaxAge = 20
#' params = log(c(120, 130, 0.1, 5))
#' Obs_delta_t=res$Obs_delta_t
#' Obs_Initlen=res$Obs_Initlen
#' Obs_Finlen=res$Obs_Finlen
#' # fit model to simulated data
#' FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
#' FittedRes$ParamEst
#' params = FittedRes$params
#' PlotOpt = 0 # all plots, 1=len inc vs init len, 2=final len vs init len, 3=resid vs init len, 4=res vs delta t,
#' # 5=annual len inc vs init len, 6=len vs age
#' PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt)
#' @export
PlotFittedTaggingGrowthModelResults <- function(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs, FittedRes, PlotOpt) {

  .pardefault <- par(no.readonly = TRUE) # store default par settings

  # if model already fitted, can input results rather than refit
  if (is.list(FittedRes)) {
    res =  FittedRes
  } else {
    res=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, MaxLen, nobs)
  }

  # params=res$params

  if (PlotOpt==0) {
    par(mfrow = c(3,2), mar=c(4,4,2,2))
  } else {
    par(mfrow = c(1,1), mar=c(4,4,2,2))
  }

  if (PlotOpt==0 | PlotOpt==1) {
    # Length inc vs init length
    xlims=Get_xaxis_scale(Obs_Initlen)
    tempdat = Obs_Finlen-Obs_Initlen
    plot(Obs_Initlen,Obs_Finlen-Obs_Initlen, cex=0.8, cex.axis=1, col="dark grey", xlab = list(" Init. Len. (mm)", cex=1.2),
         ylab = list("Len. inc. (mm)",cex=1.2), bty='n', xlim=c(0,xlims$xmax), ylim=c(0,1.2*max(tempdat)), )
    points(Obs_Initlen, res$Est_delta_Length, pch=16, col="black", cex=0.6)
  }

  # final length vs initial length
  if (PlotOpt==0 | PlotOpt==2) {
    xlims=Get_xaxis_scale(Obs_Initlen)
    xmax=xlims$xmax
    plot(Obs_Initlen,Obs_Finlen, cex=0.8, cex.axis=1, col="dark grey", xlab = list(" Init. Len. (mm)",cex=1.2),
         ylab = list("Final Len. (mm)",cex=1.2), bty='n', xlim=c(0,xlims$xmax), ylim=c(0,xlims$xmax))
    points(Obs_Initlen, res$Est_Final_Length, pch=16, col="black", cex=0.6)
  }

  # res vs initial length
  if (PlotOpt==0 | PlotOpt==3) {
    xlims=Get_xaxis_scale(Obs_Initlen)
    tempdat = Obs_Finlen-res$Est_Final_Length
    ylims=Get_yaxis_scale(tempdat)
    plot(Obs_Initlen, Obs_Finlen-res$Est_Final_Length, cex=0.8, cex.axis=1, col="dark grey", xlab = list(" Init. len (mm)",cex=1.2),
         ylab = list("Residual (mm)",cex=1.2), bty='n',xlim=c(0,xlims$xmax), ylim=c(-ylims$ymax,ylims$ymax))
    abline(h=0)
  }

  # res vs delta_t
  if (PlotOpt==0 | PlotOpt==4) {
    xlims=Get_xaxis_scale(Obs_delta_t)
    tempdat = Obs_Finlen-res$Est_Final_Length
    ylims=Get_yaxis_scale(tempdat)
    plot(Obs_delta_t, Obs_Finlen-res$Est_Final_Length, cex=0.8, cex.axis=1, col="dark grey", xlab = list(" Delta_t (days)",cex=1.2),
         ylab = list("Residual (mm)",cex=1.2), bty='n',xlim=c(0,xlims$xmax), ylim=c(-ylims$ymax,ylims$ymax))
    abline(h=0)
  }

  # estimate of annual growth increment vs specified initial length
  if (PlotOpt==0 | PlotOpt==5) {
    xlims=Get_xaxis_scale(res$plot_Init_Length)
    ylims=Get_yaxis_scale(res$plot_Length_Inc)
    plot(res$plot_Init_Length, res$plot_Length_Inc, "l", cex=0.8, cex.axis=1, xlab = list(" Init. len. (mm)",cex=1.2),
         ylab = list("Est. ann. inc. (mm)",cex=1.2), bty='n', xlim=c(0,xlims$xmax), ylim=c(0,ylims$ymax))
  }


  # estimate of length at integer ages
  if (PlotOpt==0 | PlotOpt==6) {
    xlims=Get_xaxis_scale(c(0,MaxAge))
    ylims=Get_yaxis_scale(res$plot_EstLenAtIntAge)
    plot(res$plot_IntAge, res$plot_EstLenAtIntAge, cex=0.8, cex.axis=1, pch=16, xlab = list(" Age (yrs)",cex=1.2),
         ylab = list("Length (mm)",cex=1.2), bty='n',xlim=c(0,xlims$xmax), ylim=c(0,ylims$ymax))
    lines(res$plot_IntAge,res$plot_EstLenAtIntAge, col="black")
  }

  # reset default par options
  par(.pardefault)
}


#*********************************
# Gonadosomatic index calculations
#*********************************


#' Calculate mean monthly GSIs
#'
#' This function calculates mean monthly gonaosomatic indices GSI values and associated 95 percent confidence limits
#'
#' @param MatL50 length at first maturity
#' @param FishRepdat data frame containing data for individual fish, including month of capture (Mnth), fish length (FishLen),
#' gonad wet weight (GonadWt) and fish total body weight (FishWt)
#'
#' @return data frame (Result) containing elements month (Mnth), mean GSI (Mean), standard
#' deviation (sd), monthly sample size (n), standard error (se), lower and upper 95 percent
#' confidence limits (Low95 and Up95)
#'
#' @examples
#' # Generate synthetic length composition data for a fished population at equilibrium
#' library(L3Assess)
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 15
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = NatMort
#' MaxLen = 500
#' LenInc = 1
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(150, 20) # L50, L95-L50 for gear selectivity
#' RetenParams = NA # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 350
#' vbK = 0.4
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # plot(Res$ObsDecAgeRetCatch_Fem, Res$ObsRandLenRetCatch_Fem, ylim=c(0,MaxLen), xlim=c(0,MaxAge))
#' # Now simulate length data from a Dirchlet multinomial distribution, with specified autocorrelation,
#' # to replicate non-random sampling (i.e. due to limited sampling intensity and fish schooling according to size)
#' # install.packages("dirmult")
#' library(dirmult)
#' set.seed(123)
#' nSampEvents = 5 # number of sampling events per month
#' nFishPerSampEvent = 5 # number of fish caught per sampling event
#' Mnth <- sort(rep(1:12,nSampEvents*nFishPerSampEvent))
#' nSamples = length(Mnth)
#' theta_val = 0.9 # level of autocorrelation of lengths within sampling events
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
#' FishLen <- NA
#' for (mm in 1:12) {
#'   res=SimLenFreqDat_DirMultDistn_EqMod(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#'   if (mm == 1) {
#'     FishLen <- rep(res$midpt, res$simLenFreq)
#'   } else {
#'     tempFishLen <- rep(res$midpt, res$simLenFreq)
#'     FishLen <- c(FishLen,tempFishLen)
#'   }
#' }
#' FishLen <- round(FishLen,0)
#' FishWt <- round(0.00002 * FishLen ^ 3,0)
#' # set up a linear model with specified coefficients, from which
#' # fish gonad weight data will be simulated, with error
#' a = -17.5
#' b = c(0, -0.2, -0.1, -0.1, -0.1, 0, 0.3, 0.8, 1.1, 1.5, 0.8, 0.2)
#' c = 3.5
#' ln_GonadWt <- rep(NA, nSamples)
#' for (i in 1:nSamples) {
#'   ln_GonadWt[i] =  (a + b[Mnth[i]] + c*log(FishLen[i])) + rnorm(1,0,1)
#' }
#' FishRepdat <- data.frame(GonadWt=exp(ln_GonadWt), ln_GonadWt=ln_GonadWt,
#'                          FishLen=FishLen, ln_FishLen = log(FishLen),
#'                          FishWt=FishWt, ln_FishWt = log(FishWt), Mnth)
#' # inspect simulated data
#' # par(mfrow = c(3, 2)) # Set the layout of the plots
#' # plot(FishWt ~ FishLen, data = FishRepdat, main = "Fish weight by Fish length")
#' # boxplot(FishLen ~ Mnth, data = FishRepdat, main = "FishLen by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_FishLen by Month")
#' # boxplot(GonadWt ~ Mnth, data = FishRepdat, main = "GonadWt by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_GonadWt by Month")
#' # hist(FishRepdat$FishLen)
#' # hist(FishRepdat$GonadWt)
#' # hist(FishRepdat$FishWt)
#' res=CalcMeanMonthlyGSIs(MatL50=150, FishRepdat)
#' plot(res$Mnth,res$MeanGSI, ylim=c(0,1.5*max(res$MeanGSI)))
#' points(res$Mnth,res$Low95GSI,col="blue")
#' points(res$Mnth,res$Up95GSI,col="blue")
#' @export
CalcMeanMonthlyGSIs <- function(MatL50, FishRepdat) {

  # subset data for females above at size at maturity
  subdat = FishRepdat[FishRepdat$FishLen>=MatL50,]

  # calculate GSIs for individual fish
  subdat$GSI = (subdat$GonadWt / subdat$FishWt) * 100

  # monthly means
  temp <- aggregate(GSI ~ Mnth, subdat, FUN=mean)
  Mnth = temp[,1]
  MeanGSI = temp[,2]

  # monthly counts
  temp <- aggregate(GSI ~ Mnth, subdat, FUN=length)
  nGSI = temp[,2]

  # monthly standard deviations
  temp <- aggregate(GSI ~ Mnth, subdat, FUN=sd)
  sdGSI = temp[,2]

  # create single data frame with desired fields
  results <- data.frame(Mnth=Mnth,
                       MeanGSI=MeanGSI,
                       sdGSI=sdGSI,
                       nGSI=nGSI)

  # calculate standard error
  results$seGSI <- results$sdGSI / sqrt(results$nGSI)

  # calculate lower and upper 95% CLS, assuming GSI data normally distributed
  results$Low95GSI <- results$MeanGSI - (1.96 * results$seGSI)
  results$Up95GSI <- results$MeanGSI + (1.96 * results$seGSI)

  return(results)

}

#' Calculate mean monthly, standardised gonad weights
#'
#' This function calculates mean monthly standardised gonad weights (for a specified fish length),
#' and associated 95 percent confidence limits, using GLM analysis
#'
#' @param MatL50 length at first maturity
#' @param SpecFishLength specified fish length (to which gonad weights will be standardised)
#' @param FishRepdat data frame containing data for individual fish, including month of capture (Mnth), fish length (FishLen),
#' gonad wet weight (GonadWt) and fish total body weight (FishWt)
#' @return list object containing list of calendar months with reprodutive data (MonthsWithData), monthly sample sizes,
#' output of GLM analysis, with estimated monthly effects on gonad weights (GLMResults), output of GLM prediction
#' analyses, with means and errors for standardised gonad weights (GLMPredResults), mean monthly standardised gonad weights
#' (MnthStdGonadWts_mean) and associated standard errors (MnthStdGonadWts_se), specified standard fish length inputted by user (SpecFishLength)
#' @examples
#' # Generate synthetic length composition data for a fished population at equilibrium
#' library(L3Assess)
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 15
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = NatMort
#' MaxLen = 500
#' LenInc = 1
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(150, 20) # L50, L95-L50 for gear selectivity
#' RetenParams = NA # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 350
#' vbK = 0.4
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # plot(Res$ObsDecAgeRetCatch_Fem, Res$ObsRandLenRetCatch_Fem, ylim=c(0,MaxLen), xlim=c(0,MaxAge))
#' # Now simulate length data from a Dirchlet multinomial distribution, with specified autocorrelation,
#' # to replicate non-random sampling (i.e. due to limited sampling intensity and fish schooling according to size)
#' # install.packages("dirmult")
#' library(dirmult)
#' set.seed(123)
#' nSampEvents = 5 # number of sampling events per month
#' nFishPerSampEvent = 5 # number of fish caught per sampling event
#' Mnth <- sort(rep(1:12,nSampEvents*nFishPerSampEvent))
#' nSamples = length(Mnth)
#' theta_val = 0.9 # level of autocorrelation of lengths within sampling events
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
#' FishLen <- NA
#' for (mm in 1:12) {
#'   res=SimLenFreqDat_DirMultDistn_EqMod(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#'   if (mm == 1) {
#'     FishLen <- rep(res$midpt, res$simLenFreq)
#'   } else {
#'     tempFishLen <- rep(res$midpt, res$simLenFreq)
#'     FishLen <- c(FishLen,tempFishLen)
#'   }
#' }
#' FishLen <- round(FishLen,0)
#' FishWt <- round(0.00002 * FishLen ^ 3,0)
#' # set up a linear model with specified coefficients, from which
#' # fish gonad weight data will be simulated, with error
#' a = -17.5
#' b = c(0, -0.2, -0.1, -0.1, -0.1, 0, 0.3, 0.8, 1.1, 1.5, 0.8, 0.2)
#' c = 3.5
#' ln_GonadWt <- rep(NA, nSamples)
#' for (i in 1:nSamples) {
#'   ln_GonadWt[i] =  (a + b[Mnth[i]] + c*log(FishLen[i])) + rnorm(1,0,1)
#' }
#' FishRepdat <- data.frame(GonadWt=exp(ln_GonadWt), ln_GonadWt=ln_GonadWt,
#'                          FishLen=FishLen, ln_FishLen = log(FishLen),
#'                          FishWt=FishWt, ln_FishWt = log(FishWt), Mnth)
#' # inspect simulated data
#' # par(mfrow = c(3, 2)) # Set the layout of the plots
#' # plot(FishWt ~ FishLen, data = FishRepdat, main = "Fish weight by Fish length")
#' # boxplot(FishLen ~ Mnth, data = FishRepdat, main = "FishLen by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_FishLen by Month")
#' # boxplot(GonadWt ~ Mnth, data = FishRepdat, main = "GonadWt by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_GonadWt by Month")
#' # hist(FishRepdat$FishLen)
#' # hist(FishRepdat$GonadWt)
#' # hist(FishRepdat$FishWt)
#' #
#' # plot 'true' monthly (relative) trend
#' par(mfrow=c(2,2))
#' plot(1:12, b, "o")
#' # plot GSIs
#' PlotMeanMonthlyGSIs(MatL50=150, FishRepdat, ymax=30, yint=5, GraphTitle="Females", xaxis_lab="Month",
#'                     yaxis_lab="GSI", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
#' # get monthly gonad weights, standardised for length, and plot
#' res=CalcMeanMonthlyStGonadWts(MatL50=150, SpecFishLength=300, FishRepdat)
#' res$ModelResults
#' @export
CalcMeanMonthlyStGonadWts <- function(MatL50, SpecFishLength, FishRepdat) {

  # subset data for females above at size at maturity
  subdat = FishRepdat[FishRepdat$FishLen>=MatL50,]
  MeanLenAboveL50 <- mean(subdat$FishLen)

  temp <- aggregate(FishWt  ~ Mnth, subdat, FUN=length)
  nMonthlyFishWts = temp[,2]

  # fit glm
  ln_GonadWt <- log(FishRepdat$GonadWt)
  ln_FishLen <- log(FishRepdat$FishLen)
  mod1 <- lm(ln_GonadWt ~ as.factor(Mnth) + ln_FishLen, data = subdat)
  # summary(mod1)

  # get predicted mean monthly standardised gonad weights, for fish of specified length
  newdata = data.frame(Mnth=sort(unique(subdat$Mnth)), ln_FishLen=log(SpecFishLength))
  Mod1Pred = predict(mod1, newdata, type="response", interval = "confidence")

  CalendarMonths = 1:12
  MonthsWithData = sort(unique(subdat$Mnth))
  MnthStdGonadWts_mean <- rep(NA,12)
  MnthStdGonadWts_lw95 <- rep(NA,12)
  MnthStdGonadWts_up95 <- rep(NA,12)
  MnthStdGonadWts_mean[MonthsWithData] = exp(Mod1Pred[,1])
  MnthStdGonadWts_lw95[MonthsWithData] = exp(Mod1Pred[,2])
  MnthStdGonadWts_up95[MonthsWithData] = exp(Mod1Pred[,3])

  ResultsSummary <- list(MonthsWithData = sort(unique(subdat$Mnth)),
                         MeanLenAboveL50 = MeanLenAboveL50,
                         nMonthlyFishWts = nMonthlyFishWts,
                         ModelResults = summary(mod1),
                         ModelPredictResults = Mod1Pred,
                         MnthStdGonadWts_mean = MnthStdGonadWts_mean,
                         MnthStdGonadWts_lw95 = MnthStdGonadWts_lw95,
                         MnthStdGonadWts_up95 = MnthStdGonadWts_up95,
                         SpecFishLength = SpecFishLength)

  return(ResultsSummary)

}

#' Plot mean monthly GSIs
#'
#' This function plots mean monthly GSI values and associated 95 percent confidence limits,
#' used to assess time and duration of spawning
#'
#' @param MatL50 length at first maturity
#' @param FishRepdat data frame containing data for individual fish, including month of capture (Mnth), fish length (FishLen),
#' @param ymax maximum value for y axis
#' @param yint y axis interval
#' @param GraphTitle graph title
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param SampSizelabPosAdj position of monthly sample sizes on plot above upper 95 percent CLs
#' @param SampSizelab_cex size of monthly sample size labels
#'
#' @return plot of mean month GSIs with 95 percent confidence limits
#' @examples
#' # Generate synthetic length composition data for a fished population at equilibrium
#' library(L3Assess)
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 15
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = NatMort
#' MaxLen = 500
#' LenInc = 1
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(150, 20) # L50, L95-L50 for gear selectivity
#' RetenParams = NA # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 350
#' vbK = 0.4
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # plot(Res$ObsDecAgeRetCatch_Fem, Res$ObsRandLenRetCatch_Fem, ylim=c(0,MaxLen), xlim=c(0,MaxAge))
#' # Now simulate length data from a Dirchlet multinomial distribution, with specified autocorrelation,
#' # to replicate non-random sampling (i.e. due to limited sampling intensity and fish schooling according to size)
#' # install.packages("dirmult")
#' library(dirmult)
#' set.seed(123)
#' nSampEvents = 5 # number of sampling events per month
#' nFishPerSampEvent = 5 # number of fish caught per sampling event
#' Mnth <- sort(rep(1:12,nSampEvents*nFishPerSampEvent))
#' nSamples = length(Mnth)
#' theta_val = 0.9 # level of autocorrelation of lengths within sampling events
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
#' FishLen <- NA
#' for (mm in 1:12) {
#'   res=SimLenFreqDat_DirMultDistn_EqMod(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#'   if (mm == 1) {
#'     FishLen <- rep(res$midpt, res$simLenFreq)
#'   } else {
#'     tempFishLen <- rep(res$midpt, res$simLenFreq)
#'     FishLen <- c(FishLen,tempFishLen)
#'   }
#' }
#' FishLen <- round(FishLen,0)
#' FishWt <- round(0.00002 * FishLen ^ 3,0)
#' # set up a linear model with specified coefficients, from which
#' # fish gonad weight data will be simulated, with error
#' a = -17.5
#' b = c(0, -0.2, -0.1, -0.1, -0.1, 0, 0.3, 0.8, 1.1, 1.5, 0.8, 0.2)
#' c = 3.5
#' ln_GonadWt <- rep(NA, nSamples)
#' for (i in 1:nSamples) {
#'   ln_GonadWt[i] =  (a + b[Mnth[i]] + c*log(FishLen[i])) + rnorm(1,0,1)
#' }
#' FishRepdat <- data.frame(GonadWt=exp(ln_GonadWt), ln_GonadWt=ln_GonadWt,
#'                          FishLen=FishLen, ln_FishLen = log(FishLen),
#'                          FishWt=FishWt, ln_FishWt = log(FishWt), Mnth)
#' # inspect simulated data
#' # par(mfrow = c(3, 2)) # Set the layout of the plots
#' # plot(FishWt ~ FishLen, data = FishRepdat, main = "Fish weight by Fish length")
#' # boxplot(FishLen ~ Mnth, data = FishRepdat, main = "FishLen by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_FishLen by Month")
#' # boxplot(GonadWt ~ Mnth, data = FishRepdat, main = "GonadWt by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_GonadWt by Month")
#' # hist(FishRepdat$FishLen)
#' # hist(FishRepdat$GonadWt)
#' # hist(FishRepdat$FishWt)
#' #
#' # plot 'true' monthly (relative) trend
#' par(mfrow=c(2,2))
#' plot(1:12, b, "o")
#' # plot GSIs
#' PlotMeanMonthlyGSIs(MatL50=150, FishRepdat, ymax=30, yint=5, GraphTitle="Females", xaxis_lab="Month",
#'                     yaxis_lab="GSI", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
#' @export
PlotMeanMonthlyGSIs <- function(MatL50, FishRepdat, ymax, yint, GraphTitle, xaxis_lab, yaxis_lab,
                                SampSizelabPosAdj, SampSizelab_cex) {

  res=CalcMeanMonthlyGSIs(MatL50, FishRepdat)

  if (is.na(ymax)) ymax = 2 + trunc(ceiling(max(res$Up95))/2)*2
  if (is.na(yint)) yint = 2

  MMabb = substr(month.abb, 1, 1)
  plot(res$Mnth, res$Mean, "o", pch=16, frame=F, xlim = c(1,12), ylim = c(0,ymax),
       xaxt='n', yaxt="n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab,cex=1.2), main=list(GraphTitle,cex=1.2))
  arrows(res$Mnth,res$Low95,res$Mnth,res$Up95, code=3, angle=90, length=0.03)
  AddAxesAndTickLabelsToPlot(xmin=1, xmax=12, xint=1, ymin=0, ymax, yint, cexval=1,  cexaxisval=NA, lwdval=NA,
                             lineval=0.3, lasval=2, xaxlabel = MMabb[1:12], tcklen = 0.03)
  pos = res$Up95+SampSizelabPosAdj
  text(res$Mnth, pos, res$nGSI, cex=SampSizelab_cex)
}

#' Plot mean monthly standardised gonad weights
#'
#' This function plots mean monthly standardised gonad weights (for specified fish length value) and associated
#' 95 percent confidence limits, used to assess time and duration of spawning
#'
#' @param MatL50 length at first maturity
#' @param SpecFishLength specified fish lengths to which gonad weights are to be standardised
#' @param FishRepdat data frame containing data for individual fish, including month of capture (Mnth), fish length (FishLen),
#' @param ymax maximum value for y axis
#' @param yint y axis interval
#' @param GraphTitle graph title
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param SampSizelabPosAdj position of monthly sample sizes on plot above upper 95 percent CLs
#' @param SampSizelab_cex size of monthly sample size labels
#'
#' @return plot of mean month GSIs with 95 percent confidence limits
#' @examples
#' # Generate synthetic length composition data for a fished population at equilibrium
#' library(L3Assess)
#' set.seed(123)
#' SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
#' MaxAge = 15
#' TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
#' NatMort = 4.22/MaxAge
#' FishMort = NatMort
#' MaxLen = 500
#' LenInc = 1
#' MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
#' SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
#' SelectivityVec = NA # selectivity vector
#' SelParams = c(150, 20) # L50, L95-L50 for gear selectivity
#' RetenParams = NA # L50, L95-L50 for retention
#' DiscMort = 0 # proportion of fish that die due to natural mortality
#' GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
#' Linf = 350
#' vbK = 0.4
#' CVSizeAtAge = 0.08
#' GrowthParams = c(Linf, vbK)
#' RefnceAges = NA
#' Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
#'                          SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
#' # plot(Res$ObsDecAgeRetCatch_Fem, Res$ObsRandLenRetCatch_Fem, ylim=c(0,MaxLen), xlim=c(0,MaxAge))
#' # Now simulate length data from a Dirchlet multinomial distribution, with specified autocorrelation,
#' # to replicate non-random sampling (i.e. due to limited sampling intensity and fish schooling according to size)
#' # install.packages("dirmult")
#' library(dirmult)
#' set.seed(123)
#' nSampEvents = 5 # number of sampling events per month
#' nFishPerSampEvent = 5 # number of fish caught per sampling event
#' Mnth <- sort(rep(1:12,nSampEvents*nFishPerSampEvent))
#' nSamples = length(Mnth)
#' theta_val = 0.9 # level of autocorrelation of lengths within sampling events
#' midpt = Res$midpt
#' ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
#' FishLen <- NA
#' for (mm in 1:12) {
#'   res=SimLenFreqDat_DirMultDistn_EqMod(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
#'   if (mm == 1) {
#'     FishLen <- rep(res$midpt, res$simLenFreq)
#'   } else {
#'     tempFishLen <- rep(res$midpt, res$simLenFreq)
#'     FishLen <- c(FishLen,tempFishLen)
#'   }
#' }
#' FishLen <- round(FishLen,0)
#' FishWt <- round(0.00002 * FishLen ^ 3,0)
#' # set up a linear model with specified coefficients, from which
#' # fish gonad weight data will be simulated, with error
#' a = -17.5
#' b = c(0, -0.2, -0.1, -0.1, -0.1, 0, 0.3, 0.8, 1.1, 1.5, 0.8, 0.2)
#' c = 3.5
#' ln_GonadWt <- rep(NA, nSamples)
#' for (i in 1:nSamples) {
#'   ln_GonadWt[i] =  (a + b[Mnth[i]] + c*log(FishLen[i])) + rnorm(1,0,1)
#' }
#' FishRepdat <- data.frame(GonadWt=exp(ln_GonadWt), ln_GonadWt=ln_GonadWt,
#'                          FishLen=FishLen, ln_FishLen = log(FishLen),
#'                          FishWt=FishWt, ln_FishWt = log(FishWt), Mnth)
#' # inspect simulated data
#' # par(mfrow = c(3, 2)) # Set the layout of the plots
#' # plot(FishWt ~ FishLen, data = FishRepdat, main = "Fish weight by Fish length")
#' # boxplot(FishLen ~ Mnth, data = FishRepdat, main = "FishLen by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_FishLen by Month")
#' # boxplot(GonadWt ~ Mnth, data = FishRepdat, main = "GonadWt by Month")
#' # boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_GonadWt by Month")
#' # hist(FishRepdat$FishLen)
#' # hist(FishRepdat$GonadWt)
#' # hist(FishRepdat$FishWt)
#' #
#' # plot 'true' monthly (relative) trend
#' par(mfrow=c(2,2))
#' plot(1:12, b, "o")
#' # plot GSIs
#' PlotMeanMonthlyGSIs(MatL50=150, FishRepdat, ymax=30, yint=5, GraphTitle="Females", xaxis_lab="Month",
#'                     yaxis_lab="GSI", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
#' # get monthly gonad weights, standardised for length, and plot
#' res=CalcMeanMonthlyStGonadWts(MatL50=150, SpecFishLength=300, FishRepdat)
#' res$ModelResults
#' PlotMeanMonthyStGonadWts(MatL50=150, SpecFishLength=300, FishRepdat, ymax=100, yint=20, GraphTitle="Females",
#'                          xaxis_lab="Month", yaxis_lab="Stand. gonad wt.", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
#' @export
PlotMeanMonthyStGonadWts <- function(MatL50, SpecFishLength, FishRepdat, ymax, yint, GraphTitle,
                                     xaxis_lab, yaxis_lab, SampSizelabPosAdj, SampSizelab_cex) {

  res=CalcMeanMonthlyStGonadWts(MatL50, SpecFishLength, FishRepdat)

  if (is.na(ymax)) ymax = 2 + trunc(ceiling(max(res$MnthStdGonadWts_up95))/2)*2
  if (is.na(yint)) yint = 2

  MMabb = substr(month.abb, 1, 1)
  plot(1:12, res$MnthStdGonadWts_mean, "o", pch=16, frame=F, xlim = c(1,12), ylim = c(0,ymax),
       xaxt='n', yaxt="n", xlab=list(xaxis_lab,cex=1.2), ylab=list(yaxis_lab,cex=1.2), main=list(GraphTitle,cex=1.2))
  arrows(1:12,res$MnthStdGonadWts_lw95, 1:12,res$MnthStdGonadWts_up95, code=3, angle=90, length=0.03)
  AddAxesAndTickLabelsToPlot(xmin=1, xmax=12, xint=1, ymin=0, ymax, yint, cexval=1,  cexaxisval=NA, lwdval=NA,
                             lineval=0.3, lasval=2, xaxlabel = MMabb[1:12], tcklen = 0.03)
  pos = res$MnthStdGonadWts_up95+SampSizelabPosAdj
  text(1:12, pos, res$nMonthlyFishWts, cex=SampSizelab_cex)
}

#' Calculate monthly proportions of each gonad stage
#'
#' Calculates monthly proportions of each gonad stage, used to assess time and duration of spawning
#'
#' @param MatL50 length at first maturity
#' @param FishRepdat data frame containing data for individual fish, including month of capture (MM), fish length (FishLen),

#' @return Monthly gonad stage proportions (PropGonadSt) with a list of gonad stages in data (GonadStages)
#' @examples
#' # generate some reproductive data for fish
#' set.seed(123)
#' nMnths = 12
#' MnthSampeSizes = 20
#' MMSampSize = rep(MnthSampeSizes,nMnths)
#' TotSampSize = MnthSampeSizes * nMnths
#' Mnths = 1:nMnths
#' MeanMMGSI = rep(c(rep(0.1,4),0.25, 0.7, 0.6, 0.3, rep(0.1,4)),MMSampSize)
#' randMMGSI = exp(rnorm(TotSampSize,log(MeanMMGSI),0.5))
#' MM = rep(Mnths,times=MMSampSize)
#' FishLen = round(runif(TotSampSize,100,500),0)
#' FishWt = (0.0000001 * FishLen ^ 3)
#' GonadWt = (randMMGSI / 100) * FishWt
#' FishRepdat = data.frame(MM=MM, FishLen=FishLen, GonadWt=GonadWt, FishWt=FishWt)
#' MatL50 = 200
#' JanProps = c(1,0,0,0,0,0,0,0) # specified proportions for each of the gonad stages in each month
#' FebProps = c(1,0,0,0,0,0,0,0)
#' MarProps = c(0.95,0.05,0,0,0,0,0,0)
#' AprProps = c(0.8,0.2,0,0,0,0,0,0)
#' MayProps = c(0.4,0.3,0.3,0,0,0,0,0)
#' JunProps = c(0.1,0.1,0.4,0.4,0,0,0,0)
#' JulProps = c(0.1,0,0,0.3,0.5,0,0.1,0)
#' AugProps = c(0.1,0,0,0,0,0.3,0.4,0.2)
#' SepProps = c(0.1,0,0,0,0,0,0.2,0.7)
#' OctProps = c(0.6,0,0,0,0,0,0,0.4)
#' NovProps = c(1,0,0,0,0,0,0.1,0.8)
#' DecProps = c(0.9,0,0,0,0,0,0,0.1)
#' PropsDat = data.frame(JanProps,FebProps,MarProps,AprProps,MayProps,JunProps,
#'                       JulProps,AugProps,SepProps,OctProps,NovProps,DecProps)
#' FishRepdat$GonadSt=NA
#' for (i in 1:length(FishRepdat$MM)) {
#'   mm = FishRepdat$MM[i]
#'   FishRepdat$GonadSt[i] = which(as.vector(rmultinom(1,1,PropsDat[,mm]))==1)
#' }
#' CalcMonthlyGonadStageProps(MatL50, FishRepdat)
#' @export
CalcMonthlyGonadStageProps <- function(MatL50, FishRepdat) {

  # subset data for fish above L50
  subdat = FishRepdat[FishRepdat$FishLen>=MatL50,]

  # calculate numbers of fish with each gonadal stage, caught in each month
  Table1 = table(factor(subdat$MM, levels=1:12), subdat$GonadSt)

  MinGonadStage = min(subdat$GonadSt)
  MaxGonadStage = max(subdat$GonadSt)
  GonadStages = MinGonadStage:MaxGonadStage

  # calculate monthly proportions of each gonad stage
  PropGonadSt = data.frame(matrix(nrow=12, ncol=length(GonadStages)))
  colnames(PropGonadSt) = GonadStages

  for (i in GonadStages) {
    icol = which(GonadStages==i)
    PropGonadSt[, icol] = Table1[, icol]/rowSums(Table1)
    x = which(is.nan(PropGonadSt[, icol]))
    if (length(x) > 0) {
      PropGonadSt[x,icol]=0
    }
  }

  results = list(PropGonadSt=PropGonadSt,
                 GonadStages=GonadStages)

  return(results)

}

#' Plot of calculated monthly proportions of fish gonad stages
#'
#' Plot of calculated monthly proportions of fish gonad stages, used to assess time and duration of spawning
#'
#' @param MatL50 length at first maturity
#' @param FishRepdat data frame containing data for individual fish, including month of capture (MM), and fish length (FishLen)

#' @return plot of monthly proportions of fish gonad stages
#' @examples
#' # generate some reproductive data for fish
#' set.seed(123)
#' nMnths = 12
#' MnthSampeSizes = 20
#' MMSampSize = rep(MnthSampeSizes,nMnths)
#' TotSampSize = MnthSampeSizes * nMnths
#' Mnths = 1:nMnths
#' MeanMMGSI = rep(c(rep(0.1,4),0.25, 0.7, 0.6, 0.3, rep(0.1,4)),MMSampSize)
#' randMMGSI = exp(rnorm(TotSampSize,log(MeanMMGSI),0.5))
#' MM = rep(Mnths,times=MMSampSize)
#' FishLen = round(runif(TotSampSize,100,500),0)
#' FishWt = (0.0000001 * FishLen ^ 3)
#' GonadWt = (randMMGSI / 100) * FishWt
#' FishRepdat = data.frame(MM=MM, FishLen=FishLen, GonadWt=GonadWt, FishWt=FishWt)
#' MatL50 = 200
#' JanProps = c(1,0,0,0,0,0,0,0) # specified proportions for each of the gonad stages in each month
#' FebProps = c(1,0,0,0,0,0,0,0)
#' MarProps = c(0.95,0.05,0,0,0,0,0,0)
#' AprProps = c(0.8,0.2,0,0,0,0,0,0)
#' MayProps = c(0.4,0.3,0.3,0,0,0,0,0)
#' JunProps = c(0.1,0.1,0.4,0.4,0,0,0,0)
#' JulProps = c(0.1,0,0,0.3,0.5,0,0.1,0)
#' AugProps = c(0.1,0,0,0,0,0.3,0.4,0.2)
#' SepProps = c(0.1,0,0,0,0,0,0.2,0.7)
#' OctProps = c(0.6,0,0,0,0,0,0,0.4)
#' NovProps = c(1,0,0,0,0,0,0.1,0.8)
#' DecProps = c(0.9,0,0,0,0,0,0,0.1)
#' PropsDat = data.frame(JanProps,FebProps,MarProps,AprProps,MayProps,JunProps,
#'                       JulProps,AugProps,SepProps,OctProps,NovProps,DecProps)
#' FishRepdat$GonadSt=NA
#' for (i in 1:length(FishRepdat$MM)) {
#'   mm = FishRepdat$MM[i]
#'   FishRepdat$GonadSt[i] = which(as.vector(rmultinom(1,1,PropsDat[,mm]))==1)
#' }
#' PlotMonthlyGonadStageProps(MatL50, FishRepdat)
#' @export
PlotMonthlyGonadStageProps <- function(MatL50, FishRepdat) {

  .pardefault <- par(no.readonly = TRUE) # store default par settings

  res=CalcMonthlyGonadStageProps(MatL50, FishRepdat)
  PropGonadSt=res$PropGonadSt
  GonadStages=res$GonadStages
  MMabb = substr(month.abb, 1, 1)
  par(mfcol=c(3,3), mar=c(5,4,1,1), oma=rep(0.5,4))
  for (i in GonadStages) {
    icol = which(GonadStages==i)
    plot(1:12, PropGonadSt[, icol], "o", xlim = c(1, 12), ylim = c(0, 1), cex.main = 1, main = paste("Gonad stage",i),
         bty = "n", yaxt = "n", xaxt = "n", xlab=NA, ylab=NA, lty="dotted", pch=16)
    AddAxesAndTickLabelsToPlot(xmin=1, xmax=12, xint=1, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = MMabb[1:12], tcklen = 0.03)
    if (i%in%c(1:3)) {mtext("Proportion",las=3,side=2,line=3, cex=0.8) }
    if (i%in%c(3,6,8)) {mtext("Month",las=1,side=1,line=3,adj=0.5, cex=0.8)}
  }

  # reset default par options
  par(.pardefault)
}


#*******************
# LENGTH AT MATURITY
#*******************


#' Logistic curve for length at maturity
#'
#' This function applies a logistic curve (symmetric or asymmetric) for describing probability of maturity at length
#'
#' @keywords internal
#'
#' @param params c(L50, L95) or c(Pmax, L50, L95) for CurveType=1, or c(Q, B, V) or c(Pmax, Q, B, V) for CurveType=2
#'
#' @return probability of maturity for observed lengths. Used in calculation of objective function
#' for length at maturity model
LogisticEqnLengthAtMaturity <- function(params) {

  # for estimation
  if (CurveType == 1) { # symmetric
    if (length(params)==2) { # not estimating Pmax, single sex
      L50 = params[1]; L95 = params[2]; Pmax = 1.0
    }
    if (length(params)==3) { # estimating Pmax, single sex
      Pmax = ilogit(params[1]); L50 = params[2]; L95 = params[3]
    }
    if (length(params)==4) { # not estimating Pmax, 2 sexes
      L50 = params[1:2]; L95 = params[3:4]; Pmax = c(1.0,1.0)
    }
    if (length(params)==6) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); L50 = params[3:4]; L95 = params[5:6]
    }
  } else { # asymmetric
    if (length(params)==3) { # not estimating Pmax
      Pmax = 1.0; Q=exp(params[1]); B=exp(params[2]); V=exp(params[3])
    }
    if (length(params)==4) { # estimating Pmax
      Pmax = ilogit(params[1]); Q=exp(params[2]); B=exp(params[3]); V=exp(params[4])
    }
    if (length(params)==6) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); Q = exp(params[1:2]); B = exp(params[3:4]); V = exp(params[5:6])
    }
    if (length(params)==8) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); Q = exp(params[3:4]); B = exp(params[5:6]); V = exp(params[7:8])
    }
  }

  if (nSexes==1) {
    if (CurveType == 1) { # symmetric
      results = Pmax / (1.0 + exp(- log(19) * (ObsLen - L50) / (L95 - L50)))
    } else { # asymmetric
      x = (ObsLen - B) / Q # scale length data
      results = Pmax * (1 + exp(-x)) ^ -V
    }
  }

  if (nSexes==2) {
    results = data.frame(matrix(nrow=2, ncol=length(ObsLen[1,])))
    colnames(results) = 1:length(ObsLen[1,])
    results = as.matrix(results)
    x = results
    for (i in 1:nSexes) {
      if (CurveType == 1) {
        results[i,] = Pmax[i] / (1.0 + exp(- log(19) * (ObsLen[i,] - L50[i]) / (L95[i] - L50[i])))
      } else {
        x[i,] = (ObsLen[i,] - B[i]) / Q[i] # scale length data
        results[i,] = Pmax[i] * (1 + exp(-x[i,])) ^ -V[i]
      }
    }
  }
  return(results)
}


#' Logistic curve for length at maturity (for plotting)
#'
#' This function applies a logistic curve for describing probability of maturity at length,
#' used for plotting
#'
#' @keywords internal
#'
#' @param params c(L50, L95) or c(Pmax, L50, L95) for CurveType=1, or c(Q, B, V) or c(Pmax, Q, B, V) for CurveType=2
#' @param CurveType 1 = symmetric, 2 = asymmetric
#' @param nSexes number of sexes
#' @param plotlengthrange lower and upper length class range for plotting
#'
#' @return probability of maturity for specified lengths for plotting, used for plotting maturity results
#' for length at maturity model
LogisticEqnLengthAtMaturity2 <- function(params, CurveType, nSexes, plotlengthrange) {

  # for plotting
  if (CurveType == 1) { # symmetric
    if (length(params)==2) { # not estimating Pmax, single sex
      L50 = params[1]; L95 = params[2]; Pmax = 1.0
    }
    if (length(params)==3) { # estimating Pmax, single sex
      Pmax = ilogit(params[1]); L50 = params[2]; L95 = params[3]
    }
    if (length(params)==4) { # not estimating Pmax, 2 sexes
      L50 = params[1:2]; L95 = params[3:4]; Pmax = c(1.0,1.0)
    }
    if (length(params)==6) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); L50 = params[3:4]; L95 = params[5:6]
    }
  } else { # asymmetric
    if (length(params)==3) { # not estimating Pmax
      Pmax = 1.0; Q=exp(params[1]); B=exp(params[2]); V=exp(params[3])
    }
    if (length(params)==4) { # estimating Pmax
      Pmax = ilogit(params[1]); Q=exp(params[2]); B=exp(params[3]); V=exp(params[4])
    }
    if (length(params)==6) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); Q = exp(params[1:2]); B = exp(params[3:4]); V = exp(params[5:6])
    }
    if (length(params)==8) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); Q = exp(params[3:4]); B = exp(params[5:6]); V = exp(params[7:8])
    }
  }

  if (nSexes==1) {
    plotlengths = plotlengthrange[1]:plotlengthrange[2]
    if (CurveType == 1) {
      results = Pmax / (1.0 + exp(- log(19) * (plotlengths - L50) / (L95 - L50)))
    } else {
      x = (plotlengths - B) / Q # scale length data
      results = Pmax * (1 + exp(-x)) ^ -V
    }
  }

  if (nSexes==2) {
    plotlengths_Fem = plotlengthrange[1,1]:plotlengthrange[1,2]
    plotlengths_Mal = plotlengthrange[2,1]:plotlengthrange[2,2]
    if (CurveType == 1) {
      Femresults = Pmax[1] / (1.0 + exp(- log(19) * (plotlengths_Fem - L50[1]) / (L95[1] - L50[1])))
      Malresults = Pmax[2] / (1.0 + exp(- log(19) * (plotlengths_Mal - L50[2]) / (L95[2] - L50[2])))
    } else {
      x_Fem = (plotlengths_Fem - B[1]) / Q[1] # scale length data
      Femresults = Pmax[1] * (1 + exp(-x_Fem)) ^ -V[1]
      x_Mal = (plotlengths_Mal - B[2]) / Q[2] # scale length data
      Malresults = Pmax[2] * (1 + exp(-x_Mal)) ^ -V[2]

    }

    results = list(Femresults=Femresults,
                   Malresults=Malresults)
  }

  return(results)
}

#' Logistic curve for age at maturity
#'
#' This function applies a logistic curve (symmetric or asymmetric) for describing probability of maturity at age, used
#' in calculation of objective function for maturity model
#'
#' @keywords internal
#'
#' @param params c(A50, A95) or c(Pmax, A50, A95) for CurveType=1, or c(Q, B, V) or c(Pmax, Q, B, V) for CurveType=2
#'
#' @return probability of maturity for observed ages. Used in calculation of objective function
#' for length at maturity model
LogisticEqnAgeAtMaturity <- function(params) {
  # for estimation

  if (CurveType == 1) { # symmetric
    if (length(params)==2) { # not estimating Pmax
      Pmax = 1.0; A50=params[1]; A95=params[2]
    }
    if (length(params)==3) { # estimating Pmax
      Pmax = ilogit(params[1]); A50=params[2]; A95=params[3]
    }
    if (length(params)==4) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); A50 = params[1:2]; A95 = params[3:4]
    }
    if (length(params)==6) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); A50 = params[3:4]; A95 = params[5:6]
    }
  } else { # asymmetric
    if (length(params)==3) { # not estimating Pmax
      Pmax = 1.0; Q=exp(params[1]); B=exp(params[2]); V=exp(params[3])
    }
    if (length(params)==4) { # estimating Pmax
      Pmax = ilogit(params[1]); Q=exp(params[2]); B=exp(params[3]); V=exp(params[4])
    }
    if (length(params)==6) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); Q = exp(params[1:2]); B = exp(params[3:4]); V = exp(params[5:6])
    }
    if (length(params)==8) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); Q = exp(params[3:4]); B = exp(params[5:6]); V = exp(params[7:8])
    }
  }

  if (nSexes==1) {
    if (CurveType == 1) { # symmetric
      results = Pmax / (1.0 + exp(- log(19) * (ObsAgeCl - A50) / (A95 - A50)))
    } else { # asymmetric
      x = (ObsAgeCl - B) / Q # scale age data
      results = Pmax * (1 + exp(-x)) ^ -V
    }
  }
  if (nSexes==2) {
    results = data.frame(matrix(nrow=2, ncol=length(ObsAgeCl[1,])))
    colnames(results) = 1:length(ObsAgeCl[1,])
    results = as.matrix(results)
    x = results
    for (i in 1:nSexes) {
      if (CurveType == 1) { # symmetric
        results[i,] = Pmax[i] / (1.0 + exp(-log(19) * (ObsAgeCl[i,] - A50[i]) / (A95[i] - A50[i])))
      } else { # asymmetric
        x[i,] = (ObsAgeCl[i,] - B[i]) / Q[i] # scale age data
        results[i,] = Pmax[i] * (1 + exp(-x[i,])) ^ -V[i]
      }
    }
  }
  return(results)
}


#' Logistic curve for age at maturity (for plotting)
#'
#' This function applies a logistic curve for describing probability of maturity at age, used for plotting
#' maturity results
#'
#' @keywords internal
#'
#' @param params c(A50, A95) or c(Pmax, A50, A95) for CurveType=1, or c(Q, B, V) or c(Pmax, Q, B, V) for CurveType=2
#' @param CurveType number of sexes
#' @param nSexes number of sexes
#' @param plotages age classes for plotting
#'
#' @return probability of maturity for observed ages. Used in calculation of objective function
#' for length at maturity model
LogisticEqnAgeAtMaturity2 <- function(params, CurveType, nSexes, plotages) {
  # for plotting

  if (CurveType == 1) { # symmetric
    if (length(params)==2) { # not estimating Pmax
      Pmax = 1.0; A50=params[1]; A95=params[2]
    }
    if (length(params)==3) { # estimating Pmax
      Pmax = ilogit(params[1]); A50=params[2]; A95=params[3]
    }
    if (length(params)==4) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); A50 = params[1:2]; A95 = params[3:4]
    }
    if (length(params)==6) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); A50 = params[3:4]; A95 = params[5:6]
    }
  } else { # asymmetric
    if (length(params)==3) { # not estimating Pmax
      Pmax = 1.0; Q=exp(params[1]); B=exp(params[2]); V=exp(params[3])
    }
    if (length(params)==4) { # estimating Pmax
      Pmax = ilogit(params[1]); Q=exp(params[2]); B=exp(params[3]); V=exp(params[4])
    }
    if (length(params)==6) { # not estimating Pmax, 2 sexes
      Pmax = c(1.0,1.0); Q = exp(params[1:2]); B = exp(params[3:4]); V = exp(params[5:6])
    }
    if (length(params)==8) { # estimating Pmax, 2 sexes
      Pmax = ilogit(params[1:2]); Q = exp(params[3:4]); B = exp(params[5:6]); V = exp(params[7:8])
    }
  }

  if (nSexes==1) {
    if (CurveType == 1) { # symmetric
      results = Pmax / (1.0 + exp(- log(19) * (plotages - A50) / (A95 - A50)))
    } else {
      x = (plotages - B) / Q # scale age data
      results = Pmax * (1 + exp(-x)) ^ -V
    }

  }
  if (nSexes==2) {
    results = data.frame(matrix(nrow=2, ncol=length(plotages)))
    colnames(results) = 1:length(plotages)
    results = as.matrix(results)
    x = results
    for (i in 1:nSexes) {
      if (CurveType == 1) { # symmetric
        results[i,] = Pmax[i] / (1.0 + exp(- log(19) * (plotages - A50[i]) / (A95[i] - A50[i])))
      } else {
        x[i,] = (plotages[i,] - B[i]) / Q[i] # scale age data
        results[i,] = Pmax[i] * (1 + exp(-x[i,])) ^ -V[i]
      }
    }
  }

  return(results)
}

#' Simulate length at maturity data
#'
#' Simulate length at maturity data, given specified growth, size at maturity parameters for required sample size. Can
#' simulate data for one or two sexes, and based on a symmetric or asymmetric logistic curve.
#'
#' @param nSamples number of required samples
#' @param CurveType 1 = symmetric logistic, 2 = asymmetric logistic
#' @param nSexes number of sexes
#' @param MaxAge maximum age
#' @param MinLen minimum length
#' @param MaxLen maximum length
#' @param LenInc length class interval
#' @param GrowthParams von Bertalanffy growth parameters and cv for mean size at age
#' @param MaturityParams for one or both sexes, for symmetric logistic (L50, L95, Pmax), or asymmetric logistic (Q, B, V, Pmax)
#'
#' @return Random fish lengths (ObsLen), maturity category, 1 or 0 (ObsMatCat),
#' length class lower bound (lbnd), length class mid point (midpt), proportion mature in length classes,
#' (PropMat),  plotlengths (lengths to plot maturity curve), length class sample sizes (LenClSampSize)
#'
#' @examples
#' # simulate maturity data (ignoring mortality and selectivity effects)
#' nSamples = 300
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' nSexes = 1 # single or combined sex
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 500
#' LenInc = 50
#' Linf = 300 # von Bertalanffy growth equation
#' vbK = 0.3
#' tzero = 0
#' Growth_cv = 0.08
#' GrowthParams = c(Linf, vbK, tzero, Growth_cv)
#' L50 = 180
#' L95 = 200
#' Pmax = 1.0
#' MaturityParams = c(L50,L95,Pmax)
#' Res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' # nSamples = 300
#' # CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
#' # nSexes = 1 # single or combined sex
#' # MaxAge = 20
#' # MinLen = 0
#' # MaxLen = 500
#' # LenInc = 50
#' # Linf = 300 # von Bertalanffy growth equation
#' # vbK = 0.3
#' # tzero = 0
#' # Growth_cv = 0.08
#' # GrowthParams = c(Linf, vbK, tzero, Growth_cv)
#' # Q = 20
#' # B = 200
#' # V = 2
#' # Pmax = 1.0
#' # MaturityParams = c(Q, B, V, Pmax)
#' # Res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' # nSamples = c(300, 300)
#' # nSexes = 2 # Separate sexes
#' # CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' # MaxAge = 20
#' # MinLen = 0
#' # MaxLen = 500
#' # LenInc = 50
#' # Linf = c(300, 350)
#' # vbK = c(0.3, 0.3)
#' # tzero = c(0, 0)
#' # Growth_cv = c(0.08, 0.08)
#' # GrowthParams = c(Linf, vbK, tzero, Growth_cv)
#' # L50 = c(180, 190)
#' # L95 = c(200, 210)
#' # Pmax = c(1.0, 1.0)
#' # MaturityParams = c(L50,L95,Pmax)
#' # Res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' @export
SimulateLengthAtMaturityData <- function(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams) {

  if (nSexes==1) {
    Linf = GrowthParams[1]
    vbK = GrowthParams[2]
    tzero = GrowthParams[3]
    CVSizeAtAge = GrowthParams[4]
    if (CurveType == 1) { # symmetric logistic
      L50 = MaturityParams[1]
      L95 = MaturityParams[2]
      Pmax = MaturityParams[3]
    } else { # asymmetric logistic
      Q = MaturityParams[1]
      B = MaturityParams[2]
      V = MaturityParams[3]
      Pmax = MaturityParams[4]
    }
  }
  if (nSexes==2) {
    Linf = GrowthParams[1:2]
    vbK = GrowthParams[3:4]
    tzero = GrowthParams[5:6]
    CVSizeAtAge = GrowthParams[7:8]
    if (CurveType == 1) { # symmetric logistic
      L50 = MaturityParams[1:2]
      L95 = MaturityParams[3:4]
      Pmax = MaturityParams[5:6]
    } else {
      Q = MaturityParams[1:2]
      B = MaturityParams[3:4]
      V = MaturityParams[5:6]
      Pmax = MaturityParams[7:8]
    }
  }

  # set up length bins
  lbnd = seq(MinLen,MaxLen-LenInc,LenInc)
  midpt = lbnd + LenInc/2
  nLenCl = length(lbnd)

  for (s in 1:nSexes) {
    ObsAge = runif(nSamples[s], 1, MaxAge)
    MeanLen = Linf[s] * (1 - exp(-vbK[s] * (ObsAge - tzero[s])))
    tempObsLen = rnorm(nSamples[s], MeanLen, CVSizeAtAge[s] * MeanLen)
    ObsLenCats = trunc(tempObsLen/LenInc)+1

    # assign random probabilities of maturity, given specified length
    # at maturity relationship
    tempObsLen = round(tempObsLen,0)
    nLenObs = length(tempObsLen)
    rnd = runif(nLenObs,0,1)
    tempObsMatCat = rep(NA, nLenObs)
    if (CurveType == 1) { # symmetric logistic
      ProbMat = Pmax[s] / (1 + exp(-log(19) * (tempObsLen - L50[s]) / (L95[s] - L50[s])))
    } else { # asymmetric logistic
      x = (tempObsLen - B[s]) / Q[s]
      ProbMat = Pmax[s] * (1 + exp(-x)) ^ -V[s]
    }

    for (i in 1:nLenObs) {
      if (rnd[i] <= ProbMat[i]) {
        tempObsMatCat[i] = 1
      } else {
        tempObsMatCat[i] = 0
      }
    }

    ObsFreqImm = rep(0,nLenCl)
    ObsFreqMat = rep(0,nLenCl)
    for (i in 1:nLenCl) {
      ObsFreqImm[i] <- length(which(ObsLenCats==i & tempObsMatCat==0) )
      ObsFreqMat[i] <- length(which(ObsLenCats==i & tempObsMatCat==1) )
    }
    tempLenClSampSize = ObsFreqMat + ObsFreqImm
    tempPropMat = ObsFreqMat / (ObsFreqMat + ObsFreqImm)

    # find length cut-offs for plotting
    lw = midpt[min(which(!is.nan(tempPropMat)))]
    hi = midpt[max(which(!is.nan(tempPropMat)))]
    tempplotlengths = (lw:hi)

    if (nSexes==1) {
      ObsLen = tempObsLen
      ObsMatCat = tempObsMatCat
      PropMat = tempPropMat
      plotlengths = tempplotlengths
      LenClSampSize = tempLenClSampSize
      FemObsLen = NA
      FemObsMatCat = NA
      FemPropMat = NA
      Femplotlengths = NA
      FemLenClSampSize = NA
      MalObsLen = NA
      MalObsMatCat = NA
      MalPropMat = NA
      Malplotlengths = NA
      MalLenClSampSize = NA
    }
    if (nSexes==2) {
      ObsLen = NA
      ObsMatCat = NA
      PropMat = NA
      plotlengths = NA
      LenClSampSize = NA
      if (s==1) {
        FemObsLen = tempObsLen
        FemObsMatCat = tempObsMatCat
        FemPropMat = tempPropMat
        Femplotlengths = tempplotlengths
        FemLenClSampSize = tempLenClSampSize
      } else {
        MalObsLen = tempObsLen
        MalObsMatCat = tempObsMatCat
        MalPropMat = tempPropMat
        Malplotlengths = tempplotlengths
        MalLenClSampSize = tempLenClSampSize
      }
    }
  } # sex

  results = list(ObsLen = ObsLen,
                 FemObsLen = FemObsLen,
                 MalObsLen = MalObsLen,
                 ObsMatCat = ObsMatCat,
                 FemObsMatCat = FemObsMatCat,
                 MalObsMatCat = MalObsMatCat,
                 lbnd = lbnd,
                 midpt = midpt,
                 PropMat = PropMat,
                 FemPropMat = FemPropMat,
                 MalPropMat = MalPropMat,
                 plotlengths = plotlengths,
                 Femplotlengths = Femplotlengths,
                 Malplotlengths = Malplotlengths,
                 LenClSampSize = LenClSampSize,
                 FemLenClSampSize = FemLenClSampSize,
                 MalLenClSampSize = MalLenClSampSize)

  return(results)
}

#' Simulate age at maturity data
#'
#' Simulate age at maturity data, given specified growth, size at maturity parameters and required sample size. Can simulated
#' data for 1 or 2 sexes, and based on symmetric or asymmetric logistic curve
#'
#' @param nSamples number of required samples
#' @param CurveType 1 = symmetric logistic, 2 = asymmetric logistic
#' @param nSexes number of sexes
#' @param MinAge minimum age
#' @param MaxAge maximum age
#' @param MaturityParams logistic parameters for one or two sexes, for CurveType = 1 (i.e. A50, A95) or 2 (Q, B, V, Pmax)
#'
#' @examples
#' set.seed(123)
#' MinAge = 0
#' MaxAge = 20
#' CurveType = 1  # 1 = symmetric, 2 = asymmetric
#' nSexes = 1
#' nSamples = 300
#' A50 = 4
#' A95 = 6
#' Pmax = 1.0
#' MaturityParams = c(A50, A95, Pmax)
#' SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' nSexes = 2
#' nSamples = c(300,300)
#' A50 = c(4,4.5)
#' A95 = c(6,6)
#' Pmax = c(1.0,1.0)
#' MaturityParams = c(A50, A95, Pmax)
#' SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' CurveType = 2  # 1 = symmetric, 2 = asymmetric
#' nSexes = 1
#' nSamples = 300
#' Q = 0.5
#' B = 0.5
#' V = 0.05
#' Pmax = 1.0
#' MaturityParams = c(Q,B,V,Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' # plot(res$AgeClasses, res$PropMat)
#' CurveType = 2  # 1 = symmetric, 2 = asymmetric
#' nSexes = 2
#' nSamples = c(300,300)
#' Q = c(0.5,0.5)
#' B = c(0.5,0.5)
#' V = c(0.05,0.05)
#' Pmax = c(1.0,1.0)
#' MaturityParams = c(Q,B,V,Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' # plot(res$AgeClasses, res$FemPropMat)
#' # points(res$AgeClasses, res$MalPropMat, col="blue")
#' @export
SimulateAgeAtMaturityData <- function(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams) {

  if (CurveType == 1) {  # symmetric
    if (nSexes==1) {
      A50 = MaturityParams[1]
      A95 = MaturityParams[2]
      Pmax = MaturityParams[3]
    }
    if (nSexes==2) {
      A50 = MaturityParams[1:2]
      A95 = MaturityParams[3:4]
      Pmax = MaturityParams[5:6]
    }
  }

  if (CurveType == 2) {  # asymmetric
    if (nSexes==1) {
      Q = MaturityParams[1]
      B = MaturityParams[2]
      V = MaturityParams[3]
      Pmax = MaturityParams[4]
    }
    if (nSexes==2) {
      Q = MaturityParams[1:2]
      B = MaturityParams[3:4]
      V = MaturityParams[5:6]
      Pmax = MaturityParams[7:8]
    }
  }

  AgeClasses = MinAge:MaxAge
  nAgeClasses = length(AgeClasses)

  for (s in 1:nSexes) {

    # generate random integer values between MinAge and MaxAge (ignoring selectivity, survival etc.)
    tempObsAgeCl = round(runif(nSamples[s], MinAge, MaxAge),0)
    nAgeObs = length(tempObsAgeCl)

    # assign random probabilities of maturity, given specified length
    # at maturity relationship
    rnd = runif(nAgeObs,0,1)

    tempObsMatCat = rep(NA, nAgeObs)
    if (CurveType == 1) {  # symmetric
      ProbMat = Pmax[s] / (1 + exp(-log(19) * (tempObsAgeCl - A50[s]) / (A95[s] - A50[s])))
    } else {
      x = (tempObsAgeCl - B[s]) / Q[s]
      ProbMat = Pmax[s] * (1 + exp(-x)) ^ -V[s]
    }
    for (i in 1:nAgeObs) {
      if (rnd[i] <= ProbMat[i]) {
        tempObsMatCat[i] = 1
      } else {
        tempObsMatCat[i] = 0
      }
    }

    ObsFreqImm = rep(0,nAgeClasses)
    ObsFreqMat = rep(0,nAgeClasses)
    k=0
    for (i in AgeClasses) {
      k=k+1
      ObsFreqImm[k] <- length(which(tempObsAgeCl==i & tempObsMatCat==0) )
      ObsFreqMat[k] <- length(which(tempObsAgeCl==i & tempObsMatCat==1) )
    }
    tempPropMat = ObsFreqMat / (ObsFreqMat + ObsFreqImm)
    tempAgeClSampSize = ObsFreqMat + ObsFreqImm

    # find length cut-offs for plotting
    lw = AgeClasses[min(which(!is.nan(tempPropMat)))]
    hi = AgeClasses[max(which(!is.nan(tempPropMat)))]
    tempplotages = seq(lw,hi,0.1)

    if (nSexes==1) {
      ObsAgeCl = tempObsAgeCl
      ObsMatCat = tempObsMatCat
      PropMat = tempPropMat
      plotages = tempplotages
      AgeClSampSize = tempAgeClSampSize
      FemObsAgeCl = NA
      FemObsMatCat = NA
      FemPropMat = NA
      Femplotages = NA
      FemAgeClSampSize = NA
      MalObsAgeCl = NA
      MalObsMatCat = NA
      MalPropMat = NA
      Malplotages = NA
      MalAgeClSampSize = NA
    }
    if (nSexes==2) {
      ObsAgeCl = NA
      ObsMatCat = NA
      PropMat = NA
      plotages = NA
      AgeClSampSize = NA
      if (s==1) {
        FemObsAgeCl = tempObsAgeCl
        FemObsMatCat = tempObsMatCat
        FemPropMat = tempPropMat
        Femplotages = tempplotages
        FemAgeClSampSize = tempAgeClSampSize
      }
      if (s==2) {
        MalObsAgeCl = tempObsAgeCl
        MalObsMatCat = tempObsMatCat
        MalPropMat = tempPropMat
        Malplotages = tempplotages
        MalAgeClSampSize = tempAgeClSampSize
      }
    }
  }

  results = list(ObsAgeCl = ObsAgeCl,
                 FemObsAgeCl = FemObsAgeCl,
                 MalObsAgeCl = MalObsAgeCl,
                 ObsMatCat = ObsMatCat,
                 FemObsMatCat = FemObsMatCat,
                 MalObsMatCat = MalObsMatCat,
                 AgeClasses = AgeClasses,
                 PropMat = PropMat,
                 FemPropMat = FemPropMat,
                 MalPropMat = MalPropMat,
                 plotages = plotages,
                 Femplotages = Femplotages,
                 Malplotages = Malplotages,
                 AgeClSampSize = AgeClSampSize,
                 FemAgeClSampSize = FemAgeClSampSize,
                 MalAgeClSampSize = MalAgeClSampSize)

  return(results)
}


#' Calculate negative log-likelihood for a length-based logistic maturity curve.
#'
#' Calculates the negative log-likelihood associated with a sample of fish maturity-at-length data
#' and logistic maturity parameter values.
#' @keywords internal
#'
#' @param params c(L50, L95) lengths at 50 and 95 percent maturity
#'
#' @return Negative-log likelihood associated with logistic curve fit to maturity-at-length data
CalcNLL_LogisticLengthAtMaturity <- function(params) {

  # calculate likelihood
  ProbMat = LogisticEqnLengthAtMaturity(params)

  if (nSexes==1) {
    nObs = length(ObsLen)
    Likelihood <- rep(NA,nObs)
    x=which(ObsMatCat==1)
    Likelihood[x] = ProbMat[x]
    x=which(ObsMatCat==0)
    Likelihood[x] <- 1 - ProbMat[x]
    Likelihood[which(Likelihood==0)]=1E-4
    LL <- log(Likelihood)
  }
  if (nSexes==2) {
    LL = 0
    for (i in 1:nSexes) {
      nObs = length(!is.na(ObsLen[i,]))
      Likelihood <- rep(NA,nObs)
      x=which(ObsMatCat[i,]==1)
      Likelihood[x] = ProbMat[i,x]
      x=which(ObsMatCat[i,]==0)
      Likelihood[x] <- 1 - ProbMat[i,x]
      Likelihood[which(Likelihood==0)]=1E-4
      LL = LL + log(Likelihood)
    }
  }

  # calculate the negative log-likelihood
  NLL = -sum(LL)

  # set function result to NLL
  cat("NLL",NLL,"exp(params)",exp(params),'\n')
  results = NLL

  return(results)

}

#' Calculate negative log-likelihood for an age-based logistic maturity curve.
#'
#' Calculates the negative log-likelihood associated with a sample of fish maturity-at-age data
#' and logistic maturity parameter values (or sex change data).
#' @keywords internal
#'
#' @param params c(A50, A95) ages at 50 and 95 percent maturity
#'
#' @return Negative-log likelihood associated with logistic curve fit to maturity-at-age data
CalcNLL_LogisticAgeAtMaturity <- function(params) {

  # calculate likelihood
  ProbMat = LogisticEqnAgeAtMaturity(params)

  if (nSexes==1) {
    nObs = length(ObsAgeCl)
    Likelihood <- rep(NA,nObs)
    x=which(ObsMatCat==1)
    Likelihood[x] = ProbMat[x]
    x=which(ObsMatCat==0)
    Likelihood[x] <- 1 - ProbMat[x]
    Likelihood[which(Likelihood==0)]=1E-4
    LL <- log(Likelihood)

  }
  if (nSexes==2) {
    LL = 0
    for (i in 1:nSexes) {
      nObs = length(!is.na(ObsAgeCl[i,]))
      Likelihood <- rep(NA,nObs)
      x=which(ObsMatCat[i,]==1)
      Likelihood[x] = ProbMat[i,x]
      x=which(ObsMatCat[i,]==0)
      Likelihood[x] <- 1 - ProbMat[i,x]
      Likelihood[which(Likelihood==0)]=1E-4
      LL = LL + log(Likelihood)
    }
  }

  # calculate the negative log-likelihood
  NLL = -sum(LL)

  # set function result to NLL
  cat("NLL",NLL,"params",params,'\n')
  results = NLL

  return(results)

}

#' Fit a length-based maturity curve to fish maturity-at-length data.
#'
#' This function fits a logistic maturity curve to a sample of maturity-at-length data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb
#' @keywords internal
#'
#' @param params c(L50, L95) lengths at 50 and 95 percent maturity
#' @param CurveType 1 = symmetric logistic, 2 = asymmetric logistic
#' @param nSexes number of sexes
#' @param ObsLen vector of observed lengths
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitLogisticLengthAtMaturityCurve <- function(params, CurveType, nSexes, ObsLen, ObsMatCat) {

  # run nlminb
  if (CurveType == 1) { # symmetric
    nlmb <- nlminb(params, CalcNLL_LogisticLengthAtMaturity, gradient = NULL, hessian = TRUE,
                   control=list(trace=1, eval.max=1000, iter.max=1000))
  } else { # asymmetric
    nlmb <- nlminb(params, CalcNLL_LogisticLengthAtMaturity, gradient = NULL, hessian = TRUE,
                   control=list(trace=1, eval.max=1000, iter.max=1000))
    # params=nlmb$par
    #
    # cat("fitting with Optim: Nelder-Mead",'\n')
    # nlmb <- optim(params, CalcNLL_LogisticLengthAtMaturity, method = "Nelder-Mead")
  }

  results = nlmb

  return(results)

}

#' Fit a length-based maturity curve to fish maturity-at-length data.
#'
#' This function fits a logistic maturity curve to a sample of maturity-at-age data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb
#' @keywords internal
#'
#' @param params c(A50, A95) ages at 50 and 95 percent maturity
#' @param CurveType 1 = symmetric logistic, 2 = asymmetric logistic
#' @param nSexes number of sexes
#' @param ObsAgeCl vector of observed ages
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return nlmb (stored output from internal R nlminb optimisation function)
FitLogisticAgeAtMaturityCurve <- function(params, CurveType, nSexes, ObsAgeCl, ObsMatCat) {

  # run nlminb
  if (CurveType == 1) { # symmetric
    nlmb <- nlminb(params, CalcNLL_LogisticAgeAtMaturity, gradient = NULL, hessian = TRUE,
                   control=list(trace=1, eval.max=1000, iter.max=1000))
  } else { # asymmetric
    cat("fitting with Optim: Nelder-Mead",'\n')
    nlmb <- optim(params, CalcNLL_LogisticAgeAtMaturity, method = "Nelder-Mead")
  }

  results = nlmb

  return(results)

}

#' Calculate proportions of mature fish in successive length classes
#'
#' Calculates proportions of mature fish in successive length classes,
#' and outputs some variables useful for plotting
#'
#' @param MaxLen maximum fish length
#' @param LenInc length class size increment
#' @param ObsLen vector of observed lengths
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return length class lower bounds (lbnd), length class mid points (midpt),
#' proportions of mature fish in successive length classes (PropMat), lengths
#' for plotting maturity curve (plotlengths), length class sample sizes (LenClSampSize)
#'
#' @examples
#' # generate some synthetic maturity data
#' set.seed(123)
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 400
#' LenInc = 20
#' # single sex/combined sexes
#' nSexes = 1
#' nSamples = 300
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' CVSizeAtAge = 0.1
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = 220
#' L95 = 270
#' Pmax = 1.0
#' MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' ObsLen=res$ObsLen
#' ObsMatCat=res$ObsMatCat
#' CalcPropMatureAtLength(MaxLen, LenInc, ObsLen, ObsMatCat)
#' @export
CalcPropMatureAtLength <- function(MaxLen, LenInc, ObsLen, ObsMatCat) {

  MinLen=0
  lbnd = seq(MinLen, MaxLen-LenInc, LenInc)
  midpt = lbnd + LenInc/2
  nLenCl = length(lbnd)
  ObsLenCats = trunc(ObsLen/LenInc) + 1
  ObsFreqImm = rep(0, nLenCl)
  ObsFreqMat = rep(0, nLenCl)
  for (i in 1:nLenCl) {
    ObsFreqImm[i] <- length(which(ObsLenCats == i & ObsMatCat == 0))
    ObsFreqMat[i] <- length(which(ObsLenCats == i & ObsMatCat == 1))
  }
  LenClSampSize = ObsFreqMat + ObsFreqImm
  PropMat = ObsFreqMat/(ObsFreqMat + ObsFreqImm)
  lw = midpt[min(which(!is.nan(PropMat)))]
  hi = midpt[max(which(!is.nan(PropMat)))]
  plotlengths = lw:hi

  results = list(lbnd = lbnd,
                 midpt = midpt,
                 PropMat = PropMat,
                 plotlengths = plotlengths,
                 LenClSampSize = LenClSampSize)

  return(results)
}

#' Calculate proportions of fish mature in successive age classes
#'
#' Calculates proportions of mature fish in successive age classes,
#' and outputs some variables useful for plotting
#'
#' @param MinAge minimum age class
#' @param MaxAge maximum age class
#' @param ObsAgeCl vector of observed age classes
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return Age classes (AgeClasses), proportions of mature fish in successive age classes
#' (PropMat), ages for plotting maturity curve (plotages), age class sample sizes (AgeClSampSize)
#'
#' @examples
#' # generate some synthetic maturity data
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 1
#' nSamples = 300
#' CurveType = 1
#' A50 = 4
#' A95 = 6
#' Pmax = 1.0
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' ObsAgeCl=res$ObsAgeCl
#' ObsMatCat=res$ObsMatCat
#' CalcPropMatureAtAge(MinAge, MaxAge, ObsAgeCl, ObsMatCat)
#' @export
CalcPropMatureAtAge <- function(MinAge, MaxAge, ObsAgeCl, ObsMatCat) {

  AgeClasses = MinAge:MaxAge
  nAgeClasses  = length(AgeClasses)
  ObsFreqImm = rep(0, nAgeClasses)
  ObsFreqMat = rep(0, nAgeClasses)
  k = 0
  for (i in AgeClasses) {
    k = k + 1
    ObsFreqImm[k] <- length(which(ObsAgeCl == i & ObsMatCat == 0))
    ObsFreqMat[k] <- length(which(ObsAgeCl == i & ObsMatCat == 1))
  }
  PropMat = ObsFreqMat/(ObsFreqMat + ObsFreqImm)
  AgeClSampSize = ObsFreqMat + ObsFreqImm
  lw = AgeClasses[min(which(!is.nan(PropMat)))]
  hi = AgeClasses[max(which(!is.nan(PropMat)))]
  plotages = seq(lw, hi, 0.1)
  results = list(AgeClasses = AgeClasses,
                 PropMat = PropMat,
                 plotages = plotages,
                 AgeClSampSize = AgeClSampSize)
  return(results)
}


#' Get outputs from a fitted logistic maturity curve.
#'
#' This function fits a logistic maturity curve to a sample of maturity-at-length data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' It provides various statistical outputs in include convergence statistics,
#' parameter estimated and associated 95 percent confidence limits and associated variance-covariance matrix,
#' calculated using the MASS package.
#'
#' @param params c(L50,L95) length at 50 and 95 percent maturity (L50)
#' @param nSexes c(L50,L95) number of sexes
#' @param LogisticModType 1=length-based, 2=age-based
#' @param CurveType 1=symmetric logistic, 2=asymmetric logistic
#' @param ObsLen vector of observed lengths (set to NA for age model)
#' @param ObsAgeCl vector of observed age classes (set to NA for length model)
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#' @param ErrOpt method for uncertainty calculation (1=varcov approx, 2=bootstrap)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence), sample size (SampleSize)
#' maturity parameter estimates with lower and upper 95 percent confidence limits (ParamEst), point estimates
#' for growth parameters (params), variance-covariance matrix (vcov.params)
#'
#' @examples
#' # generate synthetic length at maturity data
#' set.seed(123)
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 400
#' LenInc = 20
#' # single sex/combined sexes
#' nSexes = 1
#' nSamples = 300
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' CVSizeAtAge = 0.1
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = 220
#' L95 = 270
#' Pmax = 1.0
#' MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' # fit model to synthetic data
#' ObsAgeCl=NA
#' ObsLen=res$ObsLen
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 1 # 1=length, 2=age
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' # 2 parameter model symmetric curve
#' InitL50 = 200
#' InitL95 = 250
#' params = c(InitL50, InitL95) # without Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
#' # 3 parameter model symmetric curve
#' InitL50 = 200
#' InitL95 = 250
#' InitPmax = 0.9
#' InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
#' params = c(InitPmax_logit, InitL50, InitL95) # with Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
#' # 3 parameter model asymmetric curve
#' CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Q = 20
#' B = 200
#' V = 2
#' Pmax = 1.0
#' MaturityParams = c(Q, B, V, Pmax)
#' res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' plot(res$midpt, res$PropMat)
#' # fit model to synthetic data
#' ObsAgeCl=NA
#' ObsLen=res$ObsLen
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 1 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' InitQ = 25
#' InitB = 250
#' InitV = 1.5
#' params = c(InitQ, InitB, InitV) # with Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
#' # separate sexes
#' nSexes = 2
#' nSamples = c(500,500)
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = c(300,350)
#' vbK = c(0.3,0.3)
#' tzero = c(0,0)
#' CVSizeAtAge = c(0.1,0.1)
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = c(200,220)
#' L95 = c(240,250)
#' Pmax = c(1,1)
#' MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' ObsAgeCl=NA
#' FemObsLen=res$FemObsLen
#' MalObsLen=res$MalObsLen
#' ObsLen = as.matrix(t(data.frame(FemObsLen=FemObsLen,MalObsLen=MalObsLen)))
#' FemObsMatCat=res$FemObsMatCat
#' MalObsMatCat=res$MalObsMatCat
#' ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
#' LogisticModType = 1 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' # 2 parameter model symmetric curve
#' InitL50 = c(200, 220)
#' InitL95 = c(250, 270)
#' params = c(InitL50, InitL95) # without Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
#' # 3 parameter model symmetric curve
#' InitL50 = c(200, 220)
#' InitL95 = c(250, 270)
#' InitPmax = c(0.8, 0.8)
#' InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
#' params = c(InitPmax_logit, InitL50, InitL95) # with Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
#' # age-at-maturity
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 1
#' CurveType = 1
#' nSamples = 300
#' A50 = 4
#' A95 = 6
#' Pmax = 1.0
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' ObsLen=NA
#' ObsAgeCl=res$ObsAgeCl
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 2 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' # 2 parameter model
#' InitA50 = 3
#' InitA95 = 5
#' params = c(InitA50, InitA95) # without Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' # 3 parameter model
#' InitA50 = 4.5
#' InitA95 = 7
#' InitPmax = 0.9
#' InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
#' params = c(InitPmax_logit, InitA50, InitA95) # with Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' # # 3 parameter age-based asymmetric curve
#' # # generate synthetic length at maturity data
#' set.seed(144)
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 1
#' nSamples = 500
#' CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Q = 1 # controls spread
#' B = 5 # inflection point
#' V = 2 # controls skew or asymmetry
#' Pmax = 1
#' MaturityParams = c(Q,B,V,Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' ObsLen=NA
#' ObsAgeCl=res$ObsAgeCl
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 2 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' InitQ = 2 # related to y(0) value
#' InitB = 6 # the growth rate
#' InitV = 3 # affects near which asymptote maximum growth occurs
#' params = c(log(InitQ), InitB, log(InitV)) # without Pmax parameter
#' Res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' # separate sexes
#' nSexes = 2
#' nSamples = c(300,300)
#' CurveType = 1
#' A50 = c(4,4.5)
#' A95 = c(6,6)
#' Pmax = c(0.95,0.95)
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' FemObsAgeCl=res$FemObsAgeCl
#' MalObsAgeCl=res$MalObsAgeCl
#' ObsAgeCl = as.matrix(t(data.frame(FemObsAgeCl=FemObsAgeCl,MalObsAgeCl=MalObsAgeCl)))
#' FemObsMatCat=res$FemObsMatCat
#' MalObsMatCat=res$MalObsMatCat
#' ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
#' # 2 parameter model
#' InitA50 = c(3,3)
#' InitA95 = c(5,5)
#' params = c(InitA50, InitA95) # without Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' # 3 parameter model
#' InitA50 = c(4.5,4)
#' InitA95 = c(7,7)
#' InitPmax = c(0.9,0.9)
#' InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
#' params = c(InitPmax_logit, InitA50, InitA95) # with Pmax parameter
#' res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' @export
GetLogisticMaturityCurveResults <- function(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) {


  if (ErrOpt == 1) { # vcov approx

    # fit curve and get parameter estimates and uncertainty
    res=GetLogisticMaturityCurveResultsErrOpt1(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat)

  }

  if (ErrOpt == 2) { # bootstrap

    res=GetLogisticMaturityCurveResultsErrOpt2(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat)

  }


  results = list(nll = res$nll,
                 convergence = res$convergence,
                 SampleSize = res$SampleSize,
                 ParamEst = res$ParamEst,
                 params = res$params,
                 vcov.params = res$vcov.params,
                 ses = res$ses,
                 cor.params = res$cor.params)

  return(results)

}


#' Fit logistic maturity curve and get parameter estimates and uncertainty
#'
#' This function fits a logistic maturity curve to a sample of maturity-at-length data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' Estimates of uncertainty are derived from boostrap resampling
#'
#' @param params c(L50,L95) length at 50 and 95 percent maturity (L50)
#' @param nSexes c(L50,L95) number of sexes
#' @param LogisticModType 1=length-based, 2=age-based
#' @param CurveType 1=symmetric logistic, 2=asymmetric logistic
#' @param ObsLen vector of observed lengths (set to NA for age model)
#' @param ObsAgeCl vector of observed age classes (set to NA for length model)
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence), sample size (SampleSize)
#' maturity parameter estimates with lower and upper 95 percent confidence limits (ParamEst), point estimates
#' for growth parameters (params), variance-covariance matrix (vcov.params)

GetLogisticMaturityCurveResultsErrOpt2 <- function(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat) {

  # get maturity curve fit to initial (non-resampled) maturity data and store results
  if (LogisticModType==1) {
    nlmb = FitLogisticLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen, ObsMatCat)
  }
  if (LogisticModType==2) {
    nlmb = FitLogisticAgeAtMaturityCurve(params, CurveType, nSexes, ObsAgeCl, ObsMatCat)
  }
  nll = nlmb$objective
  convergence = nlmb$convergence
  params=nlmb$par

  # now resample (with replacement), and store parameter estimates for each resampling trial

  # Set number of bootstrap iterations
  nTrials <- 10
  ResampParams <- matrix(NA, nrow = nTrials, ncol = length(params))

  for (i in 1:nTrials) {
    cat("Trial", i, "\n")

    # Resample data
    if (LogisticModType==1) { # Length
      if (nSexes == 1) {
        FishNum_Index <- sample(1:length(ObsLen), replace = TRUE)
        ObsLen_resamp <- ObsLen[FishNum_Index]
        ObsMatCat_resamp <- ObsMatCat[FishNum_Index]
      } else if (nSexes == 2) {
        FishNum_Index <- sample(1:ncol(ObsLen), replace = TRUE)
        ObsLen_resamp <- ObsLen[, FishNum_Index]
        ObsMatCat_resamp <- ObsMatCat[, FishNum_Index]
      }
      # fit model
      fit <- tryCatch({
        nlmb = FitLogisticLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen_resamp, ObsMatCat_resamp)
      }, error = function(e) {
        cat("Error on bootstrap", i, ":", e$message, "\n")
        return(NULL)
      })
    }

    if (LogisticModType==2) { # age
      if (nSexes == 1) {
        FishNum_Index <- sample(1:length(ObsLen), replace = TRUE)
        ObsLen_resamp <- ObsLen[FishNum_Index]
        ObsMatCat_resamp <- ObsMatCat[FishNum_Index]
      } else if (nSexes == 2) {
        FishNum_Index <- sample(1:ncol(ObsLen), replace = TRUE)
        ObsLen_resamp <- ObsLen[, FishNum_Index]
        ObsMatCat_resamp <- ObsMatCat[, FishNum_Index]
      }
      # fit model
      fit <- tryCatch({
        nlmb = FitLogisticAgeAtMaturityCurve(params, CurveType, nSexes, ObsLen_resamp, ObsMatCat_resamp)
      }, error = function(e) {
        cat("Error on bootstrap", i, ":", e$message, "\n")
        return(NULL)
      })
    }

    # Save parameters if fitting was successful
    if (!is.null(fit)) {
      ResampParams[i, ] <- fit$par  # assuming your function returns a list with $par
    }
  }

  # Remove results for failed model fits
  ResampParams <- ResampParams[complete.cases(ResampParams), ]

  results = list(nll = res$nll,
                 convergence = res$convergence,
                 SampleSize = res$SampleSize,
                 ParamEst = res$ParamEst,
                 params = res$params,
                 hess.out = NA,
                 vcov.params = NA,
                 ses = NA,
                 cor.params = NA)



}




#' Fit logistic maturity curve and get parameter estimates and uncertainty
#'
#' This function fits a logistic maturity curve to a sample of maturity-at-length data
#' by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
#' Estimates of uncertainty are based on asymtotic error approximations from estimated variance
#' covariance matrix
#'
#' @param params c(L50,L95) length at 50 and 95 percent maturity (L50)
#' @param nSexes c(L50,L95) number of sexes
#' @param LogisticModType 1=length-based, 2=age-based
#' @param CurveType 1=symmetric logistic, 2=asymmetric logistic
#' @param ObsLen vector of observed lengths (set to NA for age model)
#' @param ObsAgeCl vector of observed age classes (set to NA for length model)
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#'
#' @return negative log-likelihood (nll), nlminb convergence diagnostic (convergence), sample size (SampleSize)
#' maturity parameter estimates with lower and upper 95 percent confidence limits (ParamEst), point estimates
#' for growth parameters (params), variance-covariance matrix (vcov.params)

GetLogisticMaturityCurveResultsErrOpt1 <- function(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat) {

  # fit maturity curve and get variance-covariance matrix, from fitted model, to get standard errors
  if (LogisticModType==1) {
    nlmb = FitLogisticLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen, ObsMatCat)
    hess.out = optimHess(nlmb$par, CalcNLL_LogisticLengthAtMaturity)
  }
  if (LogisticModType==2) {
    nlmb = FitLogisticAgeAtMaturityCurve(params, CurveType, nSexes, ObsAgeCl, ObsMatCat)
    hess.out = optimHess(nlmb$par, CalcNLL_LogisticAgeAtMaturity)
  }

  nll = nlmb$objective
  convergence = nlmb$convergence
  params=nlmb$par
  vcov.params = solve(hess.out)
  ses = sqrt(diag(vcov.params)) # get asymptotic standard errors of parameter estimates
  temp = diag(1/sqrt(diag(vcov.params))) # get parameter correlation matrix
  cor.params = temp %*% vcov.params %*% temp

  # calculate 95 percent confidence limits
  if (CurveType == 2) { # asymmetric logistic curve
    if (nSexes==1) {
      if (length(params)==3) { # not estimating Pmax, single sex
        EstQ = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
        EstB = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstV = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        ParamEst = t(data.frame(Q=round(EstQ,6), B=round(EstB,6), V=round(EstV,6)))
      }
      if (length(params)==4) { # estimating Pmax, single sex
        EstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
        EstQ = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        EstB = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        EstV = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        ParamEst = t(data.frame(Pmax=round(EstPmax,2), Q=round(EstQ,6), B=round(EstB,6), V=round(EstV,6)))
      }
      if (LogisticModType==1) SampleSize = length(ObsLen)
      if (LogisticModType==2) SampleSize = length(ObsAgeCl)
    } # nsexes = 1
    if (nSexes==2) {
      if (length(params)==6) { # not estimating Pmax, 2 sexes
        FemEstQ = exp(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
        MalEstQ = exp(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        MalEstB = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        FemEstB = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        MalEstV = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
        MalEstV = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))

        ParamEst = t(data.frame(FemQ=round(FemEstQ,6), FemB=round(FemEstB,6), FemV=round(FemEstV,6),
                                MalQ=round(lnMalEstQ,6), MalB=round(MalEstB,6), lnMalV=round(MalEstV,6)))
      }
      if (length(params)==6) { # estimating Pmax, 2 sexes
        FemEstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
        MalEstPmax = ilogit(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
        FemEstQ = exp(c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3]))
        MalEstQ = exp(c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4]))
        MalEstB = exp(c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5]))
        FemEstB = exp(c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6]))
        MalEstV = exp(c(nlmb$par[7], nlmb$par[7] + c(-1.96, 1.96) * ses[7]))
        MalEstV = exp(c(nlmb$par[8], nlmb$par[8] + c(-1.96, 1.96) * ses[8]))
        ParamEst = t(data.frame(FemPmax=round(FemEstPmax,2), FemQ=round(FemEstQ,6), FemB=round(FemEstB,6), FemV=round(FemEstV,6),
                                MalPmax=round(MalEstPmax,2), MalQ=round(MalEstQ,6), MalB=round(MalEstB,6), MalV=round(MalEstV,6)))
      }
      if (LogisticModType==1) SampleSize = length(which(!is.na(ObsLen[1,]))) + length(which(!is.na(ObsLen[2,])))
      if (LogisticModType==2) SampleSize = length(which(!is.na(ObsAgeCl[1,]))) + length(which(!is.na(ObsAgeCl[2,])))
    } # nSexes

  } else { # Curvetype symmetric
    if (LogisticModType==1) { # length at maturity
      if (nSexes==1) {
        if (length(params)==2) { # not estimating Pmax, single sex
          EstL50 = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
          EstL95 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          ParamEst = t(data.frame(L50=round(EstL50,2), L95=round(EstL95,2)))
        }
        if (length(params)==3) { # estimating Pmax, single sex
          EstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
          EstL50 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          EstL95 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          ParamEst = t(data.frame(Pmax=round(EstPmax,2), L50=round(EstL50,2), L95=round(EstL95,2)))
        }
        if (LogisticModType==1) SampleSize = length(ObsLen)
      } # nsexes = 1
      if (nSexes==2) {
        if (length(params)==4) { # not estimating Pmax, 2 sexes
          FemEstL50 = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
          MalEstL50 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          FemEstL95 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          MalEstL95 = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
          ParamEst = t(data.frame(FemL50=round(FemEstL50,2), FemL95=round(FemEstL95,2),
                                  MalL50=round(MalEstL50,2), MalL95=round(MalEstL95,2)))
        }
        if (length(params)==6) { # estimating Pmax, 2 sexes
          FemEstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
          MalEstPmax = ilogit(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
          FemEstL50 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          MalEstL50 = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
          FemEstL95 = c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
          MalEstL95 = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
          ParamEst = t(data.frame(FemPmax=round(FemEstPmax,2), FemL50=round(FemEstL50,2), FemL95=round(FemEstL95,2),
                                  MalPmax=round(MalEstPmax,2), MalL50=round(MalEstL50,2), MalL95=round(MalEstL95,2)))
        }
        if (LogisticModType==1) SampleSize = length(which(!is.na(ObsLen[1,]))) + length(which(!is.na(ObsLen[2,])))
      } # nsexes = 2
    } #  LogisticModType = 1
    if (LogisticModType==2) { # age at maturity
      if (nSexes==1) {
        if (length(params)==2) { # not estimating Pmax
          EstA50 = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
          EstA95 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          ParamEst = t(data.frame(A50=round(EstA50,2), A95=round(EstA95,2)))
        }
        if (length(params)==3) { # estimating Pmax
          EstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
          EstA50 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          EstA95 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          ParamEst = t(data.frame(Pmax=round(EstPmax,2), A50=round(EstA50,2), A95=round(EstA95,2)))
        }
        SampleSize = length(ObsAgeCl)
      }
      if (nSexes==2) {
        if (length(params)==4) { # not estimating Pmax, single sex
          FemEstA50 = c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1])
          MalEstA50 = c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2])
          FemEstA95 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          MalEstA95 = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
          ParamEst = t(data.frame(FemA50=round(FemEstA50,2), FemA95=round(FemEstA95,2),
                                  MalA50=round(MalEstA50,2), MalA95=round(MalEstA95,2)))
        }
        if (length(params)==6) { # estimating Pmax, 2 sexes
          FemEstPmax = ilogit(c(nlmb$par[1], nlmb$par[1] + c(-1.96, 1.96) * ses[1]))
          MalEstPmax = ilogit(c(nlmb$par[2], nlmb$par[2] + c(-1.96, 1.96) * ses[2]))
          FemEstA50 = c(nlmb$par[3], nlmb$par[3] + c(-1.96, 1.96) * ses[3])
          MalEstA50 = c(nlmb$par[4], nlmb$par[4] + c(-1.96, 1.96) * ses[4])
          FemEstA95 = c(nlmb$par[5], nlmb$par[5] + c(-1.96, 1.96) * ses[5])
          MalEstA95 = c(nlmb$par[6], nlmb$par[6] + c(-1.96, 1.96) * ses[6])
          ParamEst = t(data.frame(FemPmax=round(FemEstPmax,2), FemL50=round(FemEstA50,2), FemL95=round(FemEstA95,2),
                                  MalPmax=round(MalEstPmax,2), MalL50=round(MalEstA50,2), MalL95=round(MalEstA95,2)))
        }
        SampleSize = length(which(!is.na(ObsAgeCl[1,]))) + length(which(!is.na(ObsAgeCl[2,])))
      } # nSexes = 2
    } #  LogisticModType = 2
  } # CurveType

  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  results = list(SampleSize=SampleSize,
                 nll=nll,
                 convergence=convergence,
                 params=params,
                 vcov.params=vcov.params,
                 ses=ses,
                 cor.params=cor.params,
                 ParamEst=ParamEst)

  return(results)

}


#' Plot observed fish maturity-at-length data
#'
#' This function produces a scatter plot plot of fish length-at-maturity data
#'
#' @param nSexes number of sexes
#' @param midpt length class mid points
#' @param PropMat proportion mature in length classes
#' @param LenClSampSize length class sample sizes
#' @param xmax x axis maximum
#' @param xint x axis interval
#' @param GraphTitle graph title
#' @param SampSizelab_cex sample size graph label size
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#'
#' @return scatter plot of maturity-at-length data
#'
#' @examples
#' # generate some synthetic maturity data
#' set.seed(123)
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 400
#' LenInc = 20
#' # single sex/combined sexes
#' nSexes = 1
#' nSamples = 300
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' CVSizeAtAge = 0.1
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = 220
#' L95 = 270
#' Pmax = 1.0
#' MaturityParams = c(L50, L95, Pmax)
#' # nSexes = 2
#' # nSamples = c(300,300)
#' # CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' # Linf = c(300,350)
#' # vbK = c(0.3,0.3)
#' # tzero = c(0,0)
#' # CVSizeAtAge = c(0.1,0.1)
#' # GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' # L50 = c(220,250)
#' # L95 = c(270,290)
#' # Pmax = c(1.0,1.0)
#' # MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' midpt = res$midpt
#' # single sex
#' PropMat = res$PropMat
#' LenClSampSize=res$LenClSampSize
#' # separate sexes
#' # PropMat = t(data.frame(FemPropMat=res$FemPropMat,MalPropMat=res$MalPropMat))
#' # LenClSampSize = t(data.frame(FemLenClSampSize=res$FemLenClSampSize,MalLenClSampSize=res$MalLenClSampSize))
#' PlotLengthAtMaturityData(nSexes=1, midpt, PropMat, LenClSampSize, xmax=NA, xint=NA,
#'                          GraphTitle=NA, SampSizelab_cex=NA, xaxis_lab=NA, yaxis_lab=NA)
#' @export
PlotLengthAtMaturityData <- function(nSexes, midpt, PropMat, LenClSampSize, xmax, xint,
                                     GraphTitle, SampSizelab_cex, xaxis_lab, yaxis_lab) {

  if (is.na(xaxis_lab)) xaxis_lab = "Total length, mm"
  if (is.na(yaxis_lab)) yaxis_lab = "Prop. mature"
  if (is.na(SampSizelab_cex)) SampSizelab_cex = 0.5
  if (is.na(xmax)) xmax = 50 + trunc(ceiling(max(midpt))/50)*50
  if (is.na(xint)) xint = 50

  if (nSexes==1) {
    plot(midpt,PropMat, xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(LenClSampSize),2)
    text(midpt[temp], 1.1, LenClSampSize[temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(LenClSampSize),2)
    text(midpt[temp], 1.1, LenClSampSize[temp], cex = SampSizelab_cex, srt = 0, col="blue")
  }
  if (nSexes==2) {
    # females
    plot(midpt,PropMat[1,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(LenClSampSize[1,]),2)
    text(midpt[temp], 1.1, LenClSampSize[1,temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(LenClSampSize[1,]),2)
    text(midpt[temp], 1.1, LenClSampSize[1,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    legend("bottomright",legend="Females",bty='n', cex=0.8)
    # males
    plot(midpt,PropMat[2,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(LenClSampSize[2,]),2)
    text(midpt[temp], 1.1, LenClSampSize[2,temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(LenClSampSize[2,]),2)
    text(midpt[temp], 1.1, LenClSampSize[2,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    legend("bottomright",legend="Males",bty='n', cex=0.8)
  }
}

#' Plot observed fish maturity-at-age data
#'
#' This function produces a scatter plot plot of fish age-at-maturity data
#'
#' @param nSexes number of sexes
#' @param AgeCl vector of observed age classes
#' @param PropMat vector of calculated proportion mature by age class
#' @param AgeClSampSize vector of age class sample sizes
#' @param xmax x axis maximum
#' @param xint x axis interval
#' @param GraphTitle graph title
#' @param SampSizelab_cex sample size graph label size
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#'
#' @return scatter plot of maturity-at-age data
#'
#' @examples
#' # generate some synthetic age-at-maturity data and plot
#' # age-at-maturity data
#' MinAge = 0
#' MaxAge = 20
#' # nSexes = 1
#' # nSamples = 300
#' # A50 = 4
#' # A95 = 6
#' # Pmax = 1.0
#' # MaturityParams = c(A50, A95, Pmax)
#' nSexes = 2
#' nSamples = c(300,300)
#' CurveType = 1
#' A50 = c(4,4.5)
#' A95 = c(6,6)
#' Pmax = c(1.0,1.0)
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' AgeCl = res$AgeClasses
#' # Single sex
#' # PropMat = res$PropMat
#' # AgeClSampSize=res$AgeClSampSize
#' # separate sexes
#' PropMat = t(data.frame(FemPropMat=res$FemPropMat,MalPropMat=res$MalPropMat))
#' AgeClSampSize = t(data.frame(FemAgeClSampSize=res$FemAgeClSampSize,MalAgeClSampSize=res$MalAgeClSampSize))
#' xmax=MaxAge
#' xint=2
#' GraphTitle=NA
#' xaxis_lab=NA
#' yaxis_lab=NA
#' SampSizelab_cex=NA
#' PlotAgeAtMaturityData(nSexes, AgeCl, PropMat, AgeClSampSize, xmax, xint,
#'                       GraphTitle, SampSizelab_cex, xaxis_lab, yaxis_lab)
#' @export
PlotAgeAtMaturityData <- function(nSexes, AgeCl, PropMat, AgeClSampSize, xmax, xint,
                                  GraphTitle, SampSizelab_cex, xaxis_lab, yaxis_lab) {

  if (is.na(xmax)) xmax = 2 + trunc(ceiling(max(AgeCl))/2)*2
  if (is.na(xint)) xint = 2
  if (is.na(xaxis_lab)) xaxis_lab = "Total length, mm"
  if (is.na(yaxis_lab)) yaxis_lab = "Prop. mature"
  if (is.na(SampSizelab_cex)) SampSizelab_cex = 0.5

  if (nSexes==1) {
    plot(AgeCl,PropMat, xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(AgeClSampSize),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(AgeClSampSize),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[temp], cex = SampSizelab_cex, srt = 0, col="blue")
  }
  if (nSexes==2) {
    # females
    plot(AgeCl,PropMat[1,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(AgeClSampSize[1,]),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[1,temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(AgeClSampSize[1,]),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[1,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    legend("bottomright",legend="Females",bty='n', cex=0.8)
    # males
    plot(AgeCl,PropMat[2,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab=xaxis_lab, ylab=yaxis_lab, main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1,  cexaxisval=NA, lwdval=NA,
                               lineval=0.5, lasval=1, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    temp = seq(1,length(AgeClSampSize[2,]),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[2,temp], cex = SampSizelab_cex, srt = 0)
    temp = seq(2,length(AgeClSampSize[2,]),2)
    text(AgeCl[temp], 1.1, AgeClSampSize[2,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    legend("bottomright",legend="Males",bty='n', cex=0.8)
  }
}

#' Get confidence limits for the fitted maturity curve
#'
#' Function uses parametric resample to produce 95 percent confidence limits for a fitted maturity curve
#'
#' @param params c(L50,L95) or c(A50, A95) length or age-based logistic maturity parameters
#' @param vcov.params estimated variance-covariance matrix for maturity parameters
#' @param CurveType 1=symmetric logistic, 2=asymmetric logistic
#' @param nSexes number of sexes
#' @param LogisticModType 1=length-based, 2=age-based
#' @param plotlengthrange lower and upper length class range for plotting lengths
#' @param plotages ages for plotting if age based (set to NA is length-based)
#' @param ErrOpt 1=varcov approx, 2=boostrap
#'
#' @return scatter plot of length-at-age data
#' @export
GetConfidenceLimitsForMaturityCurve <- function(params, ErrOpt, vcov.params, CurveType, nSexes, LogisticModType, plotlengthrange, plotages) {

  # store estimated parameter distributions
  set.seed(123)
  sims = data.frame(MASS::mvrnorm(n = 500, params, vcov.params))

  if (nSexes==1) {
    if (LogisticModType==1) { # Length

      if (CurveType == 1) { # symmetric logistic
        if (length(params)==2) { names(sims) = c("L50", "L95") }
        if (length(params)==3) { names(sims) = c("Pmax", "L50", "L95") }
      } else { # asymmetric logistic
        if (length(params)==3) { names(sims) = c("Q", "B", "V") }
        if (length(params)==4) { names(sims) = c("Pmax", "Q", "B", "V") }
      }

      GetEstLenAtMat <- function(params) {
        EstLenAtMat=LogisticEqnLengthAtMaturity2(params, CurveType, nSexes, plotlengthrange)
        return(EstLenAtMat)
      }
      sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstLenAtMat)
      plotlengths = plotlengthrange[1]:plotlengthrange[2]
      sim.mat.xvals = plotlengths
    } # LogisticModType

    if (LogisticModType==2) { # Age
      if (CurveType == 1) { # symmetric logistic
        if (length(params)==2) { names(sims) = c("A50", "A95") }
        if (length(params)==3) { names(sims) = c("Pmax", "A50", "A95") }
      } else { # asymmetric logistic
        if (length(params)==3) { names(sims) = c("lnQ", "lnB", "lnV") }
        if (length(params)==4) { names(sims) = c("Pmax", "lnQ", "lnB", "lnV") }
      }

      GetEstAgeAtMat <- function(params) {
        EstAgeAtMat=LogisticEqnAgeAtMaturity2(params, CurveType, nSexes, plotages)
        return(EstAgeAtMat)
      }
      sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstAgeAtMat)
      sim.mat.xvals = plotages
    } # LogisticModType
  } # nSexes

  if (nSexes==2) { # separate sexes
    if (LogisticModType==1) { # Length
      if (CurveType == 1) { # symmetric logistic
        if (length(params)==4) { names(sims) = c("FemL50", "MalL50", "FemL95", "MalL95") }
        if (length(params)==6) { names(sims) = c("FemPmax", "MalPmax", "FemL50", "MalL50", "FemL95", "MalL95") }
      } else { # asymmetric logistic
        if (length(params)==6) { names(sims) = c("FemQ", "MalQ", "FemB", "MalB", "FemV", "MalV") }
        if (length(params)==8) { names(sims) = c("FemPmax", "MalPmax", "FemQ", "MalQ", "FemB", "MalB", "FemV", "MalV") }
      }
      GetEstLenAtMatFem <- function(params) {
        res=LogisticEqnLengthAtMaturity2(params, CurveType, nSexes, plotlengthrange)
        FemEstLenAtMat = as.vector(unlist(res$Femresults))
        return(FemEstLenAtMat)
      }
      GetEstLenAtMatMal <- function(params) {
        res=LogisticEqnLengthAtMaturity2(params, CurveType, nSexes, plotlengthrange)
        MalEstLenAtMat = as.vector(unlist(res$Malresults))
        return(MalEstLenAtMat)
      }
      Fem.sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstLenAtMatFem)
      plotlengths = plotlengthrange[1,1]:plotlengthrange[1,2]
      Fem.sim.mat.xvals = plotlengths
      Mal.sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstLenAtMatMal)
      plotlengths = plotlengthrange[2,1]:plotlengthrange[2,2]
      Mal.sim.mat.xvals = plotlengths
    }
    if (nSexes==2) {
      if (LogisticModType==2) {
        if (CurveType == 1) { # symmetric logistic
          if (length(params)==4) { names(sims) = c("FemA50", "MalA50", "FemA95", "MalA95") }
          if (length(params)==6) { names(sims) = c("FemPmax", "MalPmax", "FemA50", "MalA50", "FemA95", "MalA95") }
        } else { # asymmetric logistic
          if (length(params)==6) { names(sims) = c("FemQ", "MalQ", "FemB", "MalB", "FemV", "MalV") }
          if (length(params)==8) { names(sims) = c("FemPmax", "MalPmax", "lnFemQ", "lnMalQ", "lnFemB", "lnMalB", "lnFemV", "lnMalV") }
        }
        GetEstAgeAtMatFem <- function(params) {
          res=LogisticEqnAgeAtMaturity2(params, CurveType, nSexes, plotages)
          FemEstAgeAtMat = as.vector(unlist(res[1,]))
          return(FemEstAgeAtMat)
        }
        GetEstAgeAtMatMal <- function(params) {
          res=LogisticEqnAgeAtMaturity2(params, CurveType, nSexes, plotages)
          MalEstAgeAtMat = as.vector(unlist(res[2,]))
          return(MalEstAgeAtMat)
        }
        Fem.sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstAgeAtMatFem)
        Mal.sims.mat = apply(X=sims[,], MARGIN=1, FUN=GetEstAgeAtMatMal)
        Fem.sim.mat.xvals = plotages
        Mal.sim.mat.xvals = plotages
      }
    }
  }

  # Calculating the 2.5th an 97.5th percentile
  if (nSexes==1) {
    sim.mat.est = apply(sims.mat, 1, function(x) quantile(x, 0.5))
    sim.mat.low = apply(sims.mat, 1, function(x) quantile(x, 0.025))
    sim.mat.up = apply(sims.mat, 1, function(x) quantile(x, 0.975))
    results = list(sim.mat.est = sim.mat.est,
                   sim.mat.low = sim.mat.low,
                   sim.mat.up = sim.mat.up,
                   sim.mat.xvals = sim.mat.xvals,
                   sims.params = sims,
                   sims.curves = sims.mat)
  }
  if (nSexes==2) {
    Fem.sim.mat.est = apply(Fem.sims.mat, 1, function(x) quantile(x, 0.5))
    Fem.sim.mat.low = apply(Fem.sims.mat, 1, function(x) quantile(x, 0.025))
    Fem.sim.mat.up = apply(Fem.sims.mat, 1, function(x) quantile(x, 0.975))
    Mal.sim.mat.est = apply(Mal.sims.mat, 1, function(x) quantile(x, 0.5))
    Mal.sim.mat.low = apply(Mal.sims.mat, 1, function(x) quantile(x, 0.025))
    Mal.sim.mat.up = apply(Mal.sims.mat, 1, function(x) quantile(x, 0.975))
    results = list(Fem.sim.mat.est = Fem.sim.mat.est,
                   Fem.sim.mat.low = Fem.sim.mat.low,
                   Fem.sim.mat.up = Fem.sim.mat.up,
                   Mal.sim.mat.est = Mal.sim.mat.est,
                   Mal.sim.mat.low = Mal.sim.mat.low,
                   Mal.sim.mat.up = Mal.sim.mat.up,
                   Fem.sim.mat.xvals = Fem.sim.mat.xvals,
                   Mal.sim.mat.xvals = Mal.sim.mat.xvals,
                   sims.params = sims,
                   Fem.sims.curves = Fem.sims.mat,
                   Mal.sims.curves = Mal.sims.mat)
  }
  return(results)
}

#' Plot fitted logistic curve to fish maturity-at-length data
#'
#' This function fits a logistic maturity curve to maturity-at-length data. Various plotting
#' options provided.
#'
#' @param params observed logistic length-at-maturity parameters
#' @param CurveType 1=symmetric logistic, 2=asymmetric logistic
#' @param nSexes number of sexes
#' @param ObsLen vector of observed lengths
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#' @param LenClSampSize vector of observed length class sample sizes
#' @param midpt length class mid points
#' @param PropMat proportion mature in length classes
#' @param plotlengthrange lower and upper limits of ranges of lengths for plotting
#' @param LenClSampSize length class sample sizes
#' @param xmax x axis maximum
#' @param xint x axis interval
#' @param GraphTitle graph title
#' @param SampSizelab_cex sample size graph label size
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param SampSizelab_cex font size for sample size labels
#' @param ShowSampSizes logical (set to TRUE to show sample size labels)
#' @param PlotCLs logical (set to TRUE to show 95 percent confidence limits for curve)
#' @param ErrOpt method for uncertainty calculation (1=varcov approx, 2=bootstrap)
#'
#' @return fitted curve on scatter plot with maturity-at-length data
#' @examples
#' # plot model fitted to single sex/combined sexes
#' # generate synthetic data
#' set.seed(123)
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 400
#' LenInc = 20
#' nSexes = 1
#' nSamples = 300
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = 300
#' vbK = 0.3
#' tzero = 0
#' CVSizeAtAge = 0.1
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = 220
#' L95 = 270
#' Pmax = 1.0
#' MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' ObsAgeCl=NA
#' ObsLen=res$ObsLen
#' ObsMatCat=res$ObsMatCat
#' midpt=res$midpt
#' PropMat=res$PropMat
#' LenClSampSize=res$LenClSampSize
#' LogisticModType = 1 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' res$ParamEst
#' plotlengths=res$plotlengths
#' plotlengthrange=c(min(plotlengths),max(plotlengths))
#' # 2 parameter model
#' InitL50 = 200
#' InitL95 = 250
#' params = c(InitL50, InitL95) # without Pmax parameter
#' # # 3 parameter model
#' # InitL50 = 200
#' # InitL95 = 250
#' # InitPmax = 0.9
#' # InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
#' # params = c(InitPmax_logit, InitL50, InitL95) # with Pmax parameter
#' PlotFittedLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen, ObsMatCat, LenClSampSize, midpt, PropMat, plotlengthrange,
#'                                 xmax=400, xint=50, GraphTitle=NA, xaxis_lab="Total length, mm", yaxis_lab="Prop. mature",
#'                                 SampSizelab_cex=NA, ShowSampSizes=TRUE, ShowLegend=FALSE, PlotCLs=FALSE, ErrOpt)
#' # 3 parameter model asymmetric curve
#' CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Q = 20
#' B = 200
#' V = 2
#' Pmax = 1.0
#' MaturityParams = c(Q, B, V, Pmax)
#' res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' # fit model to synthetic data
#' ObsAgeCl=NA
#' ObsLen=res$ObsLen
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 1 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' InitQ = 25
#' InitB = 250
#' InitV = 1.5
#' params = c(InitQ, InitB, InitV) # with Pmax parameter
#' PlotFittedLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen, ObsMatCat, LenClSampSize, midpt, PropMat, plotlengthrange,
#'                                 xmax=400, xint=50, GraphTitle=NA, xaxis_lab="Total length, mm", yaxis_lab="Prop. mature",
#'                                 SampSizelab_cex=NA, ShowSampSizes=TRUE, ShowLegend=FALSE, PlotCLs=FALSE, ErrOpt)
#' # plot model fitted to separate sexes (simultaneously)
#' # generate synthetic data
#' set.seed(123)
#' MaxAge = 20
#' MinLen = 0
#' MaxLen = 400
#' LenInc = 20
#' nSexes = 2
#' nSamples = c(500,500)
#' CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Linf = c(300,350)
#' vbK = c(0.2,0.2)
#' tzero = c(0,0)
#' CVSizeAtAge = c(0.1,0.1)
#' GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
#' L50 = c(200,220)
#' L95 = c(240,250)
#' Pmax = c(1,1)
#' MaturityParams = c(L50, L95, Pmax)
#' res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
#' ObsAgeCl=NA
#' FemObsLen=res$FemObsLen
#' MalObsLen=res$MalObsLen
#' ObsLen = as.matrix(t(data.frame(FemObsLen=FemObsLen,MalObsLen=MalObsLen)))
#' FemObsMatCat=res$FemObsMatCat
#' MalObsMatCat=res$MalObsMatCat
#' ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
#' LogisticModType = 1 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' # 2 parameter model
#' InitL50 = c(200, 220)
#' InitL95 = c(250, 270)
#' params = c(InitL50, InitL95) # without Pmax parameter
#' FemPropMat=res$FemPropMat
#' MalPropMat=res$MalPropMat
#' PropMat = as.matrix(t(data.frame(FemPropMat=FemPropMat,MalPropMat=MalPropMat)))
#' FemLenClSampSize=res$FemLenClSampSize
#' MalLenClSampSize=res$MalLenClSampSize
#' LenClSampSize = as.matrix(t(data.frame(FemLenClSampSize=FemLenClSampSize,MalLenClSampSize=MalLenClSampSize)))
#' Femplotlengths = res$Femplotlengths
#' Malplotlengths = res$Malplotlengths
#' plotlengthrange = as.matrix(t(data.frame(Females=c(min(Femplotlengths),max(Femplotlengths)),
#'                                          Males=c(min(Malplotlengths),max(Malplotlengths)))))
#' midpt=res$midpt
#' PlotFittedLengthAtMaturityCurve(params, CurveType, nSexes, ObsLen, ObsMatCat, LenClSampSize, midpt, PropMat, plotlengthrange,
#'                                 xmax=400, xint=50, GraphTitle=NA, xaxis_lab="Total length, mm", yaxis_lab="Prop. mature",
#'                                 SampSizelab_cex=NA, ShowSampSizes=TRUE, ShowLegend=FALSE, PlotCLs=TRUE, ErrOpt)
#' @export
PlotFittedLengthAtMaturityCurve <- function(params, CurveType, nSexes, ObsLen, ObsMatCat, LenClSampSize, midpt, PropMat, plotlengthrange, xmax, xint,
                                            GraphTitle, xaxis_lab, yaxis_lab, SampSizelab_cex, ShowSampSizes, ShowLegend, PlotCLs, ErrOpt) {


  if (is.na(xint)) { xint = 50 }
  if (is.na(xaxis_lab)) { xaxis_lab = "Total length, mm" }
  if (is.na(yaxis_lab)) { yaxis_lab = "Prop. mature" }
  if (is.na(SampSizelab_cex)) { SampSizelab_cex = 0.5 }

  LogisticModType = 1 # length at maturity
  ObsAgeCl = NA
  res = GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt)

  if (PlotCLs == TRUE) {
    params = res$params
    vcov.params = res$vcov.params
    Res = GetConfidenceLimitsForMaturityCurve(params, ErrOpt, vcov.params, CurveType, nSexes, LogisticModType, plotlengthrange, plotages=NA)
  }

  if (nSexes==1) {
    if (is.na(xmax)) { xmax = 50 + trunc(ceiling(max(ObsLen))/50)*50 }
    plot(midpt,PropMat, xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2, las=3)
    if (ShowSampSizes) {
      temp = seq(1,length(LenClSampSize),2)
      text(midpt[temp], 1.1, LenClSampSize[temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(LenClSampSize),2)
      text(midpt[temp], 1.1, LenClSampSize[temp], cex = SampSizelab_cex, srt = 0, col="blue")
    }
    if (ShowLegend) {
      legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. L50"), lty="solid",
             bty='n', cex=0.6,lwd=-1, y.intersp=1)
    }
    plotlengths = plotlengthrange[1]:plotlengthrange[2]
    if (CurveType==1) {
      if (length(params)==2) { # not estimating Pmax
        L50 = res$ParamEst[1,1]
        L95 = res$ParamEst[2,1]
        Pmax = 1.0
      }
      if (length(params)==3) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        L50 = res$ParamEst[2,1]
        L95 = res$ParamEst[3,1]
      }
    } else {
      if (length(params)==3) { # not estimating Pmax
        Q = res$ParamEst[1,1]
        B = res$ParamEst[2,1]
        V = res$ParamEst[3,1]
        Pmax = 1.0
      }
      if (length(params)==4) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        Q = res$ParamEst[2,1]
        B = res$ParamEst[3,1]
        V = res$ParamEst[4,1]
      }
    }

    if (PlotCLs == FALSE) {
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotlengths - L50) / (L95 - L50)))
      } else {
        x = (plotlengths - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotlengths, plotprobs)
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } else { # get data for confidence limits
      lw=min(plotlengths)
      hi = max(plotlengths)
      x = c(lw:hi,hi:lw) # using shading for 95% CLs
      y = c(Res$sim.mat.low,rev(Res$sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotlengths, Res$sim.mat.est, "l", lty="solid")
      lines(plotlengths, Res$sim.mat.low, "l", lty="dotted")
      lines(plotlengths, Res$sim.mat.up, "l", lty="dotted")
      points(midpt,PropMat)
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } # else PlotCLs
  } # nSexes = 1

  if (nSexes==2) {
    # females
    if (is.na(xmax)) { xmax = 50 + trunc(ceiling(max(ObsLen[1,]))/50)*50 }
    plot(midpt,PropMat[1,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2, las=3)
    if (ShowSampSizes) {
      temp = seq(1,length(LenClSampSize[1,]),2)
      text(midpt[temp], 1.1, LenClSampSize[1,temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(LenClSampSize[1,]),2)
      text(midpt[temp], 1.1, LenClSampSize[1,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    }
    if (ShowLegend) {
      legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. L50"), lty="solid",
             bty='n', cex=0.6,lwd=-1, y.intersp=1)
      legend("bottomright",legend="Females",bty='n', cex=1)
    }

    plotlengths = plotlengthrange[1,1]:plotlengthrange[1,2]

    if (CurveType==1) {
      if (length(params)==4) { # not estimating Pmax
        L50 = res$ParamEst[1,1]
        L95 = res$ParamEst[2,1]
        Pmax = 1.0
      }
      if (length(params)==6) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        L50 = res$ParamEst[2,1]
        L95 = res$ParamEst[3,1]
      }
    } else {
      if (length(params)==3) { # not estimating Pmax
        Q = res$ParamEst[1,1]
        B = res$ParamEst[2,1]
        V = res$ParamEst[3,1]
        Pmax = 1.0
      }
      if (length(params)==4) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        Q = res$ParamEst[2,1]
        B = res$ParamEst[3,1]
        V = res$ParamEst[4,1]
      }
    }

    if (PlotCLs == FALSE) {
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotlengths - L50) / (L95 - L50)))
      } else {
        x = (plotlengths - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotlengths, plotprobs)
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } else { # get data for confidence limits
      lw=min(plotlengths)
      hi = max(plotlengths)
      x = c(lw:hi,hi:lw) # using shading for 95% CLs
      y = c(Res$Fem.sim.mat.low,rev(Res$Fem.sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotlengths, Res$Fem.sim.mat.est, "l", lty="solid")
      lines(plotlengths, Res$Fem.sim.mat.low, "l", lty="dotted")
      lines(plotlengths, Res$Fem.sim.mat.up, "l", lty="dotted")
      points(midpt,PropMat[1,])
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } # else PlotCLs

    # males
    if (is.na(xmax)) { xmax = 50 + trunc(ceiling(max(ObsLen[2,]))/50)*50 }
    plot(midpt,PropMat[2,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2, las=3)
    if (ShowSampSizes) {
      temp = seq(1,length(LenClSampSize[2,]),2)
      text(midpt[temp], 1.1, LenClSampSize[2,temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(LenClSampSize[2,]),2)
      text(midpt[temp], 1.1, LenClSampSize[2,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    }
    if (ShowLegend) {
      legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. L50"), lty="solid",
             bty='n', cex=0.6,lwd=-1, y.intersp=1)
      legend("bottomright",legend="Males",bty='n', cex=1)
    }

    plotlengths = plotlengthrange[2,1]:plotlengthrange[2,2]

    if (CurveType==1) {
      if (length(params)==4) { # not estimating Pmax
        L50 = res$ParamEst[3,1]
        L95 = res$ParamEst[4,1]
        Pmax = 1.0
      }
      if (length(params)==6) { # estimating Pmax
        Pmax = res$ParamEst[4,1]
        L50 = res$ParamEst[5,1]
        L95 = res$ParamEst[6,1]
      }
    } else {
      if (length(params)==6) { # not estimating Pmax
        Q = res$ParamEst[4,1]
        B = res$ParamEst[5,1]
        V = res$ParamEst[6,1]
        Pmax = 1.0
      }
      if (length(params)==8) { # estimating Pmax
        Pmax = res$ParamEst[5,1]
        Q = res$ParamEst[6,1]
        B = res$ParamEst[7,1]
        V = res$ParamEst[8,1]
      }
    }

    if (PlotCLs == FALSE) {
      plotlengths = plotlengthrange[2,1]:plotlengthrange[2,2]
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotlengths - L50) / (L95 - L50)))
      } else {
        x = (plotlengths - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotlengths, plotprobs)
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } else { # get data for confidence limits
      lw=min(plotlengths)
      hi = max(plotlengths)
      x = c(lw:hi,hi:lw) # using shading for 95% CLs
      y = c(Res$Mal.sim.mat.low,rev(Res$Mal.sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotlengths, Res$Mal.sim.mat.est, "l", lty="solid")
      lines(plotlengths, Res$Mal.sim.mat.low, "l", lty="dotted")
      lines(plotlengths, Res$Mal.sim.mat.up, "l", lty="dotted")
      points(midpt,PropMat[2,])
      if (CurveType==1) { points(L50,0.5*Pmax,pch=16,col="black") }
    } # else PlotCLs
  }
} # end function

#' Plot fitted logistic curve to fish maturity-at-age data.
#'
#' This function fits a logistic maturity curve to maturity-at-age data. Various plotting
#' options provided.
#'
#' @param params age-based logistic maturity curve parameters
#' @param CurveType 1=symmetric, 2=asymmetric
#' @param nSexes number of sexes
#' @param ObsAgeCl vector of observed age classes
#' @param ObsMatCat vector of observed maturity categories (0=immature, 1=mature)
#' @param AgeCl vector of observed age classes
#' @param PropMat proportion mature in age classes
#' @param AgeClSampSize vector of observed sample sizes for observed age classes
#' @param plotages ages for plotting
#' @param xmax x axis maximum
#' @param xint x axis interval
#' @param GraphTitle graph title
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param SampSizelab_cex font size for sample size labels
#' @param ShowSampSizes logical (set to TRUE to show sample size labels)
#' @param PlotCLs logical (set to TRUE to show 95 percent confidence limits for curve)
#' @param ErrOpt method for uncertainty calculation (1=varcov approx, 2=bootstrap)
#'
#' @examples
#' # plot model fitted to single sex/combined sexes
#' # generate synthetic data
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 1
#' nSamples = 500
#' CurveType = 1
#' A50 = 4
#' A95 = 6
#' Pmax = 1.0
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' AgeCl = res$AgeClasses
#' PropMat = res$PropMat
#' AgeClSampSize=res$AgeClSampSize
#' ObsAgeCl=res$ObsAgeCl
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 2 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' plotages=res$plotages
#' # 2 parameter model
#' InitA50 = 3
#' InitA95 = 5
#' params = c(InitA50, InitA95) # without Pmax parameter
#' PlotFittedAgeAtMaturityCurve(params, CurveType, nSexes, ObsAgeCl, ObsMatCat, AgeCl, PropMat, AgeClSampSize, plotages,
#'                              xmax=20, xint=2, GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, SampSizelab_cex = NA,
#'                              ShowSampSizes=FALSE, ShowLegend=FALSE, PlotCLs=FALSE, ErrOpt)
#' # plot model fitted to separate sexes (simultaneously)
#' # generate synthetic data
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 2
#' nSamples = c(500,500)
#' CurveType = 1
#' A50 = c(4,4.5)
#' A95 = c(6,6.5)
#' Pmax = c(1.0,1.0)
#' MaturityParams = c(A50, A95, Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' AgeCl = res$AgeClasses
#' FemObsAgeCl=res$FemObsAgeCl
#' MalObsAgeCl=res$MalObsAgeCl
#' ObsAgeCl = as.matrix(t(data.frame(FemObsAgeCl=FemObsAgeCl,MalObsAgeCl=MalObsAgeCl)))
#' FemObsMatCat=res$FemObsMatCat
#' MalObsMatCat=res$MalObsMatCat
#' ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
#' # 2 parameter model
#' InitA50 = c(3,3)
#' InitA95 = c(5,5)
#' params = c(InitA50, InitA95) # without Pmax parameter
#' LogisticModType=2
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' plotages=res$Femplotages # plotting routine assumes the same age range for females and males
#' FemPropMat=res$FemPropMat
#' MalPropMat=res$MalPropMat
#' PropMat = as.matrix(t(data.frame(FemPropMat=FemPropMat,MalPropMat=MalPropMat)))
#' FemAgeClSampSize=res$FemAgeClSampSize
#' MalAgeClSampSize=res$MalAgeClSampSize
#' AgeClSampSize = as.matrix(t(data.frame(FemAgeClSampSize=FemAgeClSampSize,MalAgeClSampSize=MalAgeClSampSize)))
#' PlotFittedAgeAtMaturityCurve(params, CurveType, nSexes, ObsAgeCl, ObsMatCat, AgeCl, PropMat, AgeClSampSize, plotages,
#'                              xmax=20, xint=2, GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, SampSizelab_cex = NA,
#'                              ShowSampSizes=FALSE, ShowLegend=FALSE, PlotCLs=FALSE, ErrOpt)
#' # generate synthetic age at maturity data, single sex, 3 parameter model
#' set.seed(144)
#' MinAge = 0
#' MaxAge = 20
#' nSexes = 1
#' nSamples = 500
#' CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
#' Q = 1 # controls spread
#' B = 5 # inflection point
#' V = 2 # controls skew or asymmetry
#' Pmax = 1
#' MaturityParams = c(Q,B,V,Pmax)
#' res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
#' ObsLen=NA
#' ObsAgeCl=res$ObsAgeCl
#' ObsMatCat=res$ObsMatCat
#' LogisticModType = 2 # 1=length, 2=age
#' ErrOpt = 1 # 1=varcov approx, 2=bootstrap
#' InitQ = 2 # related to y(0) value
#' InitB = 6 # the growth rate
#' InitV = 3 # affects near which asymptote maximum growth occurs
#' params = c(log(InitQ), InitB, log(InitV)) # without Pmax parameter
#' Res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
#' Res$ParamEst
#' AgeCl = res$AgeClasses
#' PropMat = res$PropMat
#' AgeClSampSize=res$AgeClSampSize
#' plotages = seq(MinAge, MaxAge, 0.1)
#' PlotFittedAgeAtMaturityCurve(params, CurveType, nSexes, ObsAgeCl, ObsMatCat, AgeCl, PropMat, AgeClSampSize, plotages,
#'                              xmax=20, xint=2, GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, SampSizelab_cex = 0.8,
#'                              ShowSampSizes=TRUE, ShowLegend=FALSE, PlotCLs=TRUE, ErrOpt)
#' @export
PlotFittedAgeAtMaturityCurve <- function(params, CurveType, nSexes, ObsAgeCl, ObsMatCat, AgeCl, PropMat, AgeClSampSize, plotages,
                                         xmax, xint, GraphTitle, xaxis_lab, yaxis_lab, SampSizelab_cex, ShowSampSizes, ShowLegend,
                                         PlotCLs, ErrOpt) {

  if (is.na(xmax)) { xmax = 2 + trunc(ceiling(max(ObsAgeCl))/2)*2 }
  if (is.na(xint)) { xint = 2 }
  if (is.na(xaxis_lab)) { xaxis_lab = "Age, y" }
  if (is.na(yaxis_lab)) { yaxis_lab = "Prop. mature" }
  if (is.na(SampSizelab_cex)) { SampSizelab_cex = 0.5 }

  # get parameter estimates
  LogisticModType = 2 # Age at maturity
  ObsLen=NA
  res = GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt)
  if (PlotCLs == TRUE) {
    # get data for confidence limits
    params = res$params
    vcov.params = res$vcov.params
    Res = GetConfidenceLimitsForMaturityCurve(params, ErrOpt, vcov.params, CurveType, nSexes, LogisticModType, plotlengthrange=NA, plotages)
  }

  # plot data
  if (nSexes==1) {
    plot(AgeCl,PropMat, xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2)
    if (ShowSampSizes) {
      temp = seq(1,length(AgeClSampSize),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(AgeClSampSize),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[temp], cex = SampSizelab_cex, srt = 0, col="blue")
      # text(AgeCl, 1.1, AgeClSampSize, cex = 0.6, srt=0)
    }
    if (ShowLegend) {
      legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. A50"), lty="solid",
             bty='n', cex=0.8,lwd=-1, y.intersp=1)
    }

    if (CurveType==1) { # symmetric logistic
      if (length(params)==2) { # not estimating Pmax
        A50 = res$ParamEst[1,1]
        A95 = res$ParamEst[2,1]
        Pmax = 1.0
      }
      if (length(params)==3) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        A50 = res$ParamEst[2,1]
        A95 = res$ParamEst[3,1]
      }
    } else {
      if (length(params)==3) { # not estimating Pmax
        Q = res$ParamEst[1,1]
        B = res$ParamEst[2,1]
        V = res$ParamEst[3,1]
        Pmax = 1.0
      }
      if (length(params)==4) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        Q = res$ParamEst[2,1]
        B = res$ParamEst[3,1]
        V = res$ParamEst[4,1]
      }
    }

    if (PlotCLs == FALSE) {
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotages - A50) / (A95 - A50)))
      } else {
        x = (plotages - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotages, plotprobs)
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } else {
      # plot confidence limits
      x1 = plotages
      x2 = rev(plotages)
      x = c(x1,x2) # using shading for 95% CLs
      y = c(Res$sim.mat.low,rev(Res$sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotages, Res$sim.mat.est, "l", lty="solid")
      lines(plotages, Res$sim.mat.low, "l", lty="dotted")
      lines(plotages, Res$sim.mat.up, "l", lty="dotted")
      points(AgeCl,PropMat)
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } # else PlotCLs
  }
  if (nSexes==2) {

    # females
    plot(AgeCl,PropMat[1,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2)
    if (ShowSampSizes) {
      temp = seq(1,length(AgeClSampSize[1,]),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[1,temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(AgeClSampSize[1,]),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[1,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    }
    if (ShowLegend) {
      legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. L50"), lty="solid",
             bty='n', cex=0.6,lwd=-1, y.intersp=1)
      legend("bottomright",legend="Females",bty='n', cex=1)
    }

    if (CurveType==1) { # symmetric logistic
      if (length(params)==4) { # not estimating Pmax
        A50 = res$ParamEst[1,1]
        A95 = res$ParamEst[2,1]
        Pmax = 1.0
      }
      if (length(params)==6) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        A50 = res$ParamEst[2,1]
        A95 = res$ParamEst[3,1]
      }
    } else {
      if (length(params)==6) { # not estimating Pmax
        Q = res$ParamEst[1,1]
        B = res$ParamEst[2,1]
        V = res$ParamEst[3,1]
        Pmax = 1.0
      }
      if (length(params)==8) { # estimating Pmax
        Pmax = res$ParamEst[1,1]
        Q = res$ParamEst[2,1]
        B = res$ParamEst[3,1]
        V = res$ParamEst[4,1]
      }
    }

    if (PlotCLs == FALSE) {
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotages - A50) / (A95 - A50)))
      } else {
        x = (plotages - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotages, plotprobs)
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } else {
      # plot confidence limits
      x1 = plotages
      x2 = rev(plotages)
      x = c(x1,x2) # using shading for 95% CLs
      y = c(Res$Fem.sim.mat.low,rev(Res$Fem.sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotages, Res$Fem.sim.mat.est, "l", lty="solid")
      lines(plotages, Res$Fem.sim.mat.low, "l", lty="dotted")
      lines(plotages, Res$Fem.sim.mat.up, "l", lty="dotted")
      points(AgeCl,PropMat[1,])
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } # else PlotCLs

    # males
    plot(AgeCl,PropMat[2,], xlim=c(0,xmax), ylim=c(0,1.1),
         frame=F, xaxt='n', yaxt='n', xlab="", ylab="", main=GraphTitle)
    AddAxesAndTickLabelsToPlot(xmin=0, xmax, xint, ymin=0, ymax=1, yint=0.2, cexval=1.2,  cexaxisval=1, lwdval=NA,
                               lineval=0, lasval=3, xaxlabel = seq(0,xmax,xint), tcklen = 0.03)
    mtext(xaxis_lab, 1, line = 2.5, cex=1.2)
    mtext(yaxis_lab, 2, line = 2.5, cex=1.2)
    if (ShowSampSizes) {
      temp = seq(1,length(AgeClSampSize[2,]),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[2,temp], cex = SampSizelab_cex, srt = 0)
      temp = seq(2,length(AgeClSampSize[2,]),2)
      text(AgeCl[temp], 1.1, AgeClSampSize[2,temp], cex = SampSizelab_cex, srt = 0, col="blue")
    }

    if (ShowLegend) {
    legend(0.01*xmax,0.975, pch=c(1,16), legend=c("Observed", "Est. L50"), lty="solid",
           bty='n', cex=0.6,lwd=-1, y.intersp=1)
    legend("bottomright",legend="Males",bty='n', cex=1)
    }

    if (CurveType==1) { # symmetric logistic
      if (length(params)==4) { # not estimating Pmax
        A50 = res$ParamEst[3,1]
        A95 = res$ParamEst[4,1]
        Pmax = 1.0
      }
      if (length(params)==6) { # estimating Pmax
        Pmax = res$ParamEst[4,1]
        A50 = res$ParamEst[5,1]
        A95 = res$ParamEst[6,1]
      }
    } else {
      if (length(params)==6) { # not estimating Pmax
        Q = res$ParamEst[4,1]
        B = res$ParamEst[5,1]
        V = res$ParamEst[6,1]
        Pmax = 1.0
      }
      if (length(params)==8) { # estimating Pmax
        Pmax = res$ParamEst[5,1]
        Q = res$ParamEst[6,1]
        B = res$ParamEst[7,1]
        V = res$ParamEst[8,1]
      }
    }

    if (PlotCLs == FALSE) {
      if (CurveType==1) {
        plotprobs = Pmax / (1.0 + exp(- log(19) * (plotages - A50) / (A95 - A50)))
      } else {
        x = (plotages - B) / Q # scale age data
        plotprobs = Pmax * (1 + exp(-x)) ^ -V
      }
      lines(plotages, plotprobs)
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } else {
      # plot confidence limits
      x1 = plotages
      x2 = rev(plotages)
      x = c(x1,x2) # using shading for 95% CLs
      y = c(Res$Mal.sim.mat.low,rev(Res$Mal.sim.mat.up))
      polygon(x,y,col="light grey",border=NA)
      lines(plotages, Res$Mal.sim.mat.est, "l", lty="solid")
      lines(plotages, Res$Mal.sim.mat.low, "l", lty="dotted")
      lines(plotages, Res$Mal.sim.mat.up, "l", lty="dotted")
      points(AgeCl,PropMat[2,])
      if (CurveType==1) { points(A50,0.5*Pmax,pch=16,col="black")  }
    } # else PlotCLs
  }
} # end function




# ***************************
# Weight-length relationships
# ***************************

#' Simulate weight-length data
#'
#' Function simulates weight-length data for individual fish
#'
#' @param nsamples number of required samples
#' @param minlen minimum length
#' @param maxlen maximum length
#' @param lenwt_a weight length equation parameter
#' @param lenwt_b weight length equation parameter
#' @param lenwt_sd standard deviation of data (in log space)
#'
#' @return randomly-generated fish lengths (FishLen) and fish weights (FishWt)
#'
#' @examples
#' # generate synthetic data
#' set.seed(123)
#' nsamples = 100
#' minlen = 20
#' maxlen = 500
#' lenwt_a = 0.00002
#' lenwt_b = 3
#' lenwt_sd = 0.1
#' Res = SimulateWeightLengthData(nsamples, minlen, maxlen,
#'                                lenwt_a, lenwt_b, lenwt_sd)
#' Res
#' @export
SimulateWeightLengthData <- function(nsamples, minlen, maxlen,
                                     lenwt_a, lenwt_b, lenwt_sd) {

  plotlen = seq(minlen, maxlen, 1)
  plotwt = lenwt_a * plotlen ^ lenwt_b
  rand_len = round(runif(nsamples, minlen, maxlen),0)
  exp_wt = lenwt_a * rand_len ^ lenwt_b
  rand_err = rnorm(nsamples, 0, lenwt_sd)
  rand_wt = exp_wt * exp(rand_err) * exp(0.5 * lenwt_sd ^ 2)

  results = list(FishLen=rand_len,
                 FishWt=rand_wt)
  return(results)

}

#' Fit weight-length relationship and get parameter estimates
#'
#' Function fits a weight-length relationship to data in log space using the lm() functin in R
#'
#' @param FishLen vector of observed fish lengths
#' @param FishWt vector of observed fish weights
#'
#' @return outputs of linear model used to fit weight-length relationship in log space (mod1), weight-length parameters
#' with associated 95 percent confidence limits (ParamEst), adjusted coefficient of determination for fitted model, R2
#' (adj.r.squared), lengths for plotting line in log space (plot_lnlen), weights for plotting line in log space (plot_lnwt),
#' 95 percent confidence limits for fitted line (conf_int), 95 percent prediction intervals (pred_int), residual variance (Resid_var)
#'
#' @examples
#' # generate synthetic data
#' nsamples = 100
#' minlen = 20
#' maxlen = 500
#' lenwt_a = 0.00002
#' lenwt_b = 3
#' lenwt_sd = 0.1
#' Res = SimulateWeightLengthData(nsamples, minlen, maxlen,
#'                                lenwt_a, lenwt_b, lenwt_sd)
#' FishLen=Res$FishLen
#' FishWt=Res$FishWt
#' # Fit weight-length relationship and get results
#' GetWeightLengthRegressionResults(FishLen, FishWt)
#' @export
GetWeightLengthRegressionResults <- function(FishLen, FishWt) {

  # do linear regression, using lm function
  ln_len = log(FishLen)
  ln_wt = log(FishWt)
  dat = data.frame(ln_len, ln_wt)
  mod1 = lm(ln_wt ~ ln_len, data=dat)
  out=summary(mod1)

  # get parameter estimates and 95% confidence limits
  log_a = out$coefficients[1,1]
  log_a_se = out$coefficients[1,2]
  b = out$coefficients[2,1]
  b_se = out$coefficients[2,2]
  ln_lenwt_a = c(log_a, log_a + c(-1.96,1.96) * log_a_se)
  lenwt_b = c(b, b + c(-1.96,1.96) * b_se)
  ParamEst = t(data.frame(ln_lenwt_a=round(ln_lenwt_a,3), lenwt_b=round(lenwt_b,3)))
  colnames(ParamEst) = c("Estimate","lw_95%CL","up_95%CL")

  # get R2
  adj.r.squared = out$adj.r.squared

  # get sample size
  sampsize = length(FishLen)

  # get residual variance
  Resid_var = (summary(mod1)$sigma)^2

  # get estimated weights, in log space, for sequential (logged) fish length values within range of x data
  max_lnlen = max(ln_len)
  min_lnlen = min(ln_len)
  plot_lnlen = seq(min_lnlen, max_lnlen, 0.1)
  plot_lnwt = (log_a + b * plot_lnlen)

  # get 95% confidence limits
  newdat <- data.frame(ln_len=plot_lnlen)
  conf_int = predict(mod1, newdata = newdat, interval = 'confidence')
  pred_int = predict(mod1, newdata = newdat, interval = 'predict')

  results = list(mod1 = mod1,
                 ParamEst = ParamEst,
                 adj.r.squared = adj.r.squared,
                 plot_lnlen = plot_lnlen,
                 plot_lnwt = plot_lnwt,
                 conf_int = conf_int,
                 pred_int = pred_int,
                 Resid_var=Resid_var)

  return(results)
}

#' Plot weight length relationship over data in normal space
#'
#' @param FishLen vector of fish lengths
#' @param FishWt vector of fish weights
#' @param xmin minimum value for x axis
#' @param ymin minimum value for y axis
#' @param xmax maximum value for x axis
#' @param ymax maximum value for y axis
#' @param xint interval for x axis
#' @param yint interval for y axis
#' @param GraphTitle title for graph
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param PlotCLs logical (set to TRUE to plot 95 percent confidence and prediction limits)
#'
#' @return fitted curve on scatter plot with weight-length data
#' @examples
#' # generate synthetic data
#' nsamples = 100
#' minlen = 20
#' maxlen = 500
#' lenwt_a = 0.00002
#' lenwt_b = 3
#' lenwt_sd = 0.1
#' Res = SimulateWeightLengthData(nsamples, minlen, maxlen,
#'                                lenwt_a, lenwt_b, lenwt_sd)
#' FishLen=Res$FishLen
#' FishWt=Res$FishWt/1000 # kg
#' # Fit weight-length relationship and plot results (in normal space)
#' PlotWeightLengthRel_NormalSpace(FishLen, FishWt, xmin=NA, ymin=NA, xmax=NA, ymax=NA, xint=NA, yint=NA,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T)
#' @export
PlotWeightLengthRel_NormalSpace <- function(FishLen, FishWt, xmin, ymin, xmax, ymax, xint, yint,
                                            GraphTitle, xaxis_lab, yaxis_lab, PlotCLs) {

  # get regression analysis results
  res = GetWeightLengthRegressionResults(FishLen, FishWt)

  # get default axis limits and intervals
  ylims=Get_yaxis_scale(FishWt)
  xlims=Get_xaxis_scale(FishLen)

  if (is.na(xmin)) xmin = 0
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymin)) ymin = 0
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "Weight (kg)"

  # plot in normal space (with bias correction applied)
  plot(FishLen, FishWt, pch=16, cex=0.5, bty='n', ylim=c(ymin,ymax), xlim=c(xmin,xmax), xaxt='n', yaxt='n',
       ylab="", xlab="", main=GraphTitle, cex.main=1)
  AddAxesAndTickLabelsToPlot(xmin, xmax, xint, ymin, ymax, yint, cexval=1.2,  cexaxisval=1, lwdval=1.75,
                             lineval=0, lasval=1, xaxlabel = seq(xmin,xmax,xint), tcklen = 0.03)
  mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2,lwd=1.75)
  mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2,lwd=1.75)
  lines(exp(res$plot_lnlen), exp(res$conf_int[,1]) * exp(0.5*res$Resid_var), col="black", lty="solid")

  if (PlotCLs == T) {
    lines(exp(res$plot_lnlen), exp(res$conf_int[,2]) * exp(0.5*res$Resid_var), col="dark grey", lty="dashed")
    lines(exp(res$plot_lnlen), exp(res$conf_int[,3]) * exp(0.5*res$Resid_var), col="dark grey", lty="dashed")
    lines(exp(res$plot_lnlen), exp(res$pred_int[,2]) * exp(0.5*res$Resid_var), col="blue", lty="dotted")
    lines(exp(res$plot_lnlen), exp(res$pred_int[,3]) * exp(0.5*res$Resid_var), col="blue", lty="dotted")
    legend(0.05*xmax, 0.95*ymax, pch=c(1,16), legend=c("Estimate", "low 95% CL", "up 95% CL", "low 95% PL", "up 95% PL"),
           lty=c("solid","dashed","dashed","dotted","dotted"),
           col=c("black","dark grey","dark grey","blue","blue"),
           bty='n', cex=0.8,lwd=1, y.intersp=1)
  }
}

#' Plot weight length relationship over data in log space
#'
#' @param FishLen vector of fish lengths
#' @param FishWt vector of fish weights
#' @param xmin minimum value for x axis
#' @param ymin minimum value for y axis
#' @param xmax maximum value for x axis
#' @param ymax maximum value for y axis
#' @param xint interval for x axis
#' @param yint interval for y axis
#' @param GraphTitle title for graph
#' @param xaxis_lab x axis label
#' @param yaxis_lab y axis label
#' @param PlotCLs logical (set to TRUE to plot 95 percent confidence and prediction limits)
#'
#' @return fitted curve on scatter plot with weight-length data in log space
#' @examples
#' # generate synthetic data
#' nsamples = 100
#' minlen = 20
#' maxlen = 500
#' lenwt_a = 0.00002
#' lenwt_b = 3
#' lenwt_sd = 0.1
#' Res = SimulateWeightLengthData(nsamples, minlen, maxlen,
#'                                lenwt_a, lenwt_b, lenwt_sd)
#' FishLen=Res$FishLen
#' FishWt=Res$FishWt/1000 # kg
#' # Fit weight-length relationship and plot results (in log space)
#' PlotWeightLengthRel_LogSpace(FishLen, FishWt, xmin=NA, ymin=NA, xmax=NA, ymax=NA, xint=NA, yint=NA,
#' GraphTitle=NA, xaxis_lab=NA, yaxis_lab=NA, PlotCLs=T)
#' @export
PlotWeightLengthRel_LogSpace <- function(FishLen, FishWt, xmin, ymin, xmax, ymax, xint, yint,
                                         GraphTitle, xaxis_lab, yaxis_lab, PlotCLs) {

  # get regression analysis results
  res = GetWeightLengthRegressionResults(FishLen, FishWt)

  # get default axis limits and intervals
  ylims=Get_yaxis_scale(log(FishWt))
  xlims=Get_xaxis_scale(log(FishLen))

  if (is.na(xmin)) xmin = xlims$xmin
  if (is.na(xmax)) xmax = xlims$xmax
  if (is.na(xint)) xint = xlims$xint
  if (is.na(ymin)) ymin = ylims$ymin
  if (is.na(ymax)) ymax = ylims$ymax
  if (is.na(yint)) yint = ylims$yint
  if (is.na(xaxis_lab)) xaxis_lab = "ln Length (mm)"
  if (is.na(yaxis_lab)) yaxis_lab = "ln Weight (kg)"

  # plot in log space
  plot(log(FishLen), log(FishWt), pch=16, cex=0.5, bty='n', ylim=c(ymin,ymax), xlim=c(xmin,xmax), xaxt='n', yaxt='n',
       ylab="", xlab="", main=GraphTitle, cex.main=1)
  AddAxesAndTickLabelsToPlot(xmin, xmax, xint, ymin, ymax, yint, cexval=1.2,  cexaxisval=1, lwdval=1.75,
                             lineval=0, lasval=1, xaxlabel = seq(xmin,xmax,xint), tcklen = 0.03)
  mtext(yaxis_lab,las=3,side=2,line=3,cex=1.2,lwd=1.75)
  mtext(xaxis_lab,las=1,side=1,line=3,cex=1.2,lwd=1.75)
  lines(res$plot_lnlen, res$conf_int[,1], col="black", lty="solid")

  if (PlotCLs == T) {
    lines(res$plot_lnlen, res$conf_int[,2], col="dark grey", lty="dashed")
    lines(res$plot_lnlen, res$conf_int[,3], col="dark grey", lty="dashed")
    lines(res$plot_lnlen, res$pred_int[,2], col="blue", lty="dotted")
    lines(res$plot_lnlen, res$pred_int[,3], col="blue", lty="dotted")
    legend(1.05*xmin, 0.95*ymax, pch=c(1,16), legend=c("Estimate", "low 95% CL", "up 95% CL", "low 95% PL", "up 95% PL"),
           lty=c("solid","dashed","dashed","dotted","dotted"),
           col=c("black","dark grey","dark grey","blue","blue"),
           bty='n', cex=0.8,lwd=1, y.intersp=1)
  }
}


#' Simulate multiple years of length at age data with time-varying growth
#'
#' This function simulates multiple connsecutive years of length-at-age data
#' with time-varying growth, based on year specific von Bertalanffy growth
#' parameters (Linf and/or vbK). Data are simulated assuming constant mortality
#' and selectivity.
#'
#' @param Linf asymptotic length (by sex, and year, if multiple years)
#' @param vbK growth coefficient (by sex, and year, if multiple years)
#' @param tzero combined sex, single value
#' @param CV coefficient of variation, constant for mean length at age
#' @param nLinfVals number of consecutive years of Linf params
#' @param nvbKVals number of consecutive years of vbK params
#' @param StartYear first year of length-at-age data
#' @param EndYear last year of length-at-age data
#' @param minIntAge minimum integer age
#' @param MaxIntAge maximum integer age
#' @param nSexes number of sexes
#' @return results (RandLenAtAge, EstLenAtAge, EstLenAtIntAge, nYrs, nIntAges, nObs, LenAtAgeDat)
#' @examples
#' library(L3Assess)
#' library(RTMB)
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort = 0.11
#' FMort = 0.11
#' ZMort = FMort + NatMort
#' SelA50 = c(7,6)
#' SelA95 = c(8,7)
#' SampleSize = 300 # required sample size. For 2 sex model, same sample size generated for each sex.
#' Res=SimAgeFreqData_EqMod(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = Res$CatchSample
#' par(mfrow=c(1,1),mar=c(4,3,2,2))
#' plot(MinAge:MaxAge, ObsAgeFreq[1,],col="red",'o')
#' lines(MinAge:MaxAge, ObsAgeFreq[2,],col="blue")
#' # Specify information required to generate length-at-age data with
#' # time-varying growth
#' StartYear = 2000
#' EndYear = 2020
#' nYrs = EndYear - StartYear + 1
#' minIntAge = 0
#' MaxIntAge = 40
#' nSexes = 2
#' # Specify true growth parameter values. Linf and/or vbK can be year-specific.
#' nLinfVals = 1
#' nvbKVals = nYrs
#' Linf_Fem = rep(1000,nLinfVals)
#' Linf_Mal = rep(1100,nLinfVals)
#' Linf = data.frame(data.frame(Linf_Fem=Linf_Fem,Linf_Mal=Linf_Mal))
#' vbK_Fem = seq(0.1,0.2,0.1/nYrs)
#' vbK_Mal = seq(0.1,0.2,0.1/nYrs)
#' vbK = data.frame(data.frame(vbK_Fem=vbK_Fem,vbK_Mal=vbK_Mal))
#' tzero = 0 # Specify common tzero for both sexes to simplify
#' CV = 0.1
#' Res = SimulateTimeVaryingLengthAtAgeData(Linf, vbK, tzero, CV, nLinfVals, nvbKVals, StartYear, EndYear, minIntAge, MaxIntAge, nSexes)
#' @export
SimulateTimeVaryingLengthAtAgeData <- function(Linf, vbK, tzero, CV, nLinfVals, nvbKVals, StartYear, EndYear, minIntAge, MaxIntAge, nSexes) {

  nYrs = EndYear - StartYear + 1
  nIntAges = length(minIntAge:MaxIntAge)
  nObs = SampleSize * nYrs * nSexes

  if (nLinfVals==1) {
    Linf_Fem = rep(Linf_Fem,nYrs)
    Linf_Mal = rep(Linf_Mal,nYrs)
    Linf = data.frame(Linf_Fem=Linf_Fem,Linf_Mal=Linf_Mal)
  }
  if (nvbKVals==1) {
    vbK_Fem = rep(vbK_Fem,nYrs)
    vbK_Mal = rep(vbK_Mal,nYrs)
    vbK = data.frame(vbK_Fem=vbK_Fem, vbK_Mal=vbK_Mal)
  }


  # generate data with with year, sex, integer age and decimal age
  YrIndex = sort(rep(1:nYrs,SampleSize*nSexes))
  Sex = rep(sort(rep(1:2,SampleSize)),nYrs)
  IntAge = rep(c(rep(MinAge:MaxAge, ObsAgeFreq[1,]),rep(MinAge:MaxAge, ObsAgeFreq[2,])),nYrs)
  Age = IntAge + runif(nObs,0,1)
  ObsLenAtAgeDat <- data.frame(YrIndex=YrIndex, Sex=Sex, Age=Age)

  Res = CalcYearEffectGrowthDynamics(ObsLenAtAgeDat, Linf, vbK, tzero)

  nObs = length(ObsLenAtAgeDat$Age)
  RandLenAtAge = rep(0,nObs)
  for (j in 1:nObs) {
    SD = CV * Res$EstLenAtAge[j]
    RandLenAtAge[j] = round(Res$EstLenAtAge[j] + rnorm(1,0,SD),0)
  }

  ObsLenAtAgeDat$RandLenAtAge = RandLenAtAge

  results = list(RandLenAtAge = RandLenAtAge,
                 EstLenAtAge = Res$EstLenAtAge,
                 EstLenAtIntAge = Res$EstLenAtIntAge,
                 nYrs=nYrs,
                 nIntAges=nIntAges,
                 nObs=nObs,
                 ObsLenAtAgeDat=ObsLenAtAgeDat)

  return(results)
}

#' Simulate multiple years of length at age data with time-varying growth
#'
#' This function simulates multiple connsecutive years of length-at-age data
#' with time-varying growth, based on year specific von Bertalanffy growth
#' parameters (Linf and/or vbK). Data are simulated assuming constant mortality
#' and selectivity.
#'
#' @keywords internal
#'
#' @param ObsLenAtAgeDat year specific length-at-age data
#' @param Linf asymptotic length (by sex, and year, if multiple years)
#' @param vbK growth coefficient (by sex, and year, if multiple years)
#' @param tzero combined sex, single value
#'
#' @return results (EstLenAtIntAge, EstLenAtAge)
CalcYearEffectGrowthDynamics <- function(ObsLenAtAgeDat, Linf, vbK, tzero) {

  nSexes = length(unique(ObsLenAtAgeDat$Sex))
  MinYear = min(ObsLenAtAgeDat$YrIndex)
  MaxYear = max(ObsLenAtAgeDat$YrIndex)
  nYrs = MaxYear - MinYear + 1
  IntAge = floor(ObsLenAtAgeDat$Age)
  MinAge = 0
  MaxAge = max(IntAge) + 1
  nIntAges = MaxAge - MinAge
  nObs = length(ObsLenAtAgeDat$Age)

  # used to prevent possibility of negative growth
  eps = 0.1
  slope = log(19)/eps
  bound = 0.001
  EstLenAtIntAge <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  tempAnnGrowthInc <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  EstAnnGrowthInc <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  for (y in 1:nYrs) {  # Year
    for (s in 1:nSexes) { # Sex
      for (j in 1:nIntAges) {  # Integer ages
        i = minIntAge + j - 1  # age
        if (y == 1) {
          if (j==1) {
            EstAnnGrowthInc[y,s,j] = 0
            EstLenAtIntAge[y,s,j] = Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))
          } else {
            EstAnnGrowthInc[y,s,j] = (Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))) - EstLenAtIntAge[y,s,j-1]
            EstLenAtIntAge[y,s,j] = EstLenAtIntAge[y,s,j-1] + EstAnnGrowthInc[y,s,j]
          }
        }
        if (y > 1) {
          if (j == 1) {
            EstAnnGrowthInc[y,s,j] = 0
            EstLenAtIntAge[y,s,j] = Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))
          }
          if (j > 1) {
            # prevent possibility of negative growth
            tempAnnGrowthInc[y,s,j] = ((Linf[y,s] - EstLenAtIntAge[y-1,s,j-1]) * (1 - exp(-vbK[y,s])))
            EstAnnGrowthInc[y,s,j] = 1/(1+exp(slope*(bound-tempAnnGrowthInc[y,s,j]))) * tempAnnGrowthInc[y,s,j]
            EstLenAtIntAge[y,s,j] = EstLenAtIntAge[y-1,s,j-1] + EstAnnGrowthInc[y,s,j]
          } # j>1
        } # y>1
      } # j
    } # s
  } # y

  EstLenAtAge = rep(0,nObs)
  tempGrowthInc = rep(0,nObs)
  EstGrowthInc = rep(0,nObs)
  for (j in 1:nObs) {

    y = ObsLenAtAgeDat$YrIndex[j]
    s = ObsLenAtAgeDat$Sex[j]
    i = floor(ObsLenAtAgeDat$Age[j])
    a = ObsLenAtAgeDat$Age[j]

    tempGrowthInc[j] = (Linf[y,s] - EstLenAtIntAge[y,s,i]) * (1 - exp(-vbK[y,s] * (a - i)))
    EstGrowthInc[j] = 1/(1+exp(slope*(bound-tempGrowthInc[j]))) * tempGrowthInc[j]
    EstLenAtAge[j] = EstLenAtIntAge[y,s,i] + EstGrowthInc[j]

  } # j

  results = list(EstLenAtIntAge = EstLenAtIntAge,
                 EstLenAtAge = EstLenAtAge,
                 EstGrowthInc = EstGrowthInc)

  return(results)

}

#' Calculate NLL associated with fit of year-effects growth model
#'
#' This function calculates the negative log-likelihood associated with the fit
#' of a von Bertalanffy growth model with year-effects parameters (Linf and/or vbK).
#' The function links to RTMB.
#'
#' @keywords internal
#'
#' @param params von Bertalanffy growth parameters (some of which are time-varying)
#'
#' @return NLL (and several Report and ADreport variables, associated with RTMB)
Calculate_NLL_YearEffectGrowthMod <- function(params) {

  # RTMB
  # https://groups.google.com/g/tmb-users/c/HlPqkfcCa1g
  "[<-" <- ADoverload("[<-")

  # get parameters
  RTMB::getAll(dat, params, warn=FALSE)
  Linf_Fem <- exp(lnLinf_Fem)
  Linf_Mal <- exp(lnLinf_Mal)
  vbK_Fem <- exp(lnvbK_Fem)
  vbK_Mal <- exp(lnvbK_Mal)

  if (nLinfVals==1) {
    Linf_Fem = rep(Linf_Fem,nYrs)
    Linf_Mal = rep(Linf_Mal,nYrs)
  }
  if (nvbKVals==1) {
    vbK_Fem = rep(vbK_Fem,nYrs)
    vbK_Mal = rep(vbK_Mal,nYrs)
  }
  Linf = data.frame(Linf_Fem=Linf_Fem,Linf_Mal=Linf_Mal)
  vbK = data.frame(vbK_Fem=vbK_Fem, vbK_Mal=vbK_Mal)

  # get estimated lengths at specified integer ages
  nSexes = length(unique(ObsLenAtAgeDat$Sex))
  MinYear = min(ObsLenAtAgeDat$YrIndex)
  MaxYear = max(ObsLenAtAgeDat$YrIndex)
  nYrs = MaxYear - MinYear + 1
  IntAge = floor(ObsLenAtAgeDat$Age)
  MinAge = 0
  MaxAge = max(IntAge) + 1
  nIntAges = MaxAge - MinAge
  nObs = length(ObsLenAtAgeDat$Age)
  # used to prevent possibility of negative growth
  eps = 0.1
  slope = log(19)/eps
  bound = 0.001
  EstLenAtIntAge <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  tempAnnGrowthInc <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  EstAnnGrowthInc <- array(0, dim=(c(nYrs, nSexes, nIntAges)))
  for (y in 1:nYrs) {  # Year
    for (s in 1:nSexes) { # Sex
      for (j in 1:nIntAges) {  # Integer ages
        i = minIntAge + j - 1  # age
        if (y == 1) {
          if (j==1) {
            EstAnnGrowthInc[y,s,j] = 0
            EstLenAtIntAge[y,s,j] = Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))
          } else {
            EstAnnGrowthInc[y,s,j] = (Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))) - EstLenAtIntAge[y,s,j-1]
            EstLenAtIntAge[y,s,j] = EstLenAtIntAge[y,s,j-1] + EstAnnGrowthInc[y,s,j]
          }
        }
        if (y > 1) {
          if (j == 1) {
            EstAnnGrowthInc[y,s,j] = 0
            EstLenAtIntAge[y,s,j] = Linf[y,s] * (1 - exp(-vbK[y,s] * (i - tzero)))
          }
          if (j > 1) {
            # prevent possibility of negative growth
            tempAnnGrowthInc[y,s,j] = ((Linf[y,s] - EstLenAtIntAge[y-1,s,j-1]) * (1 - exp(-vbK[y,s])))
            EstAnnGrowthInc[y,s,j] = 1/(1+exp(slope*(bound-tempAnnGrowthInc[y,s,j]))) * tempAnnGrowthInc[y,s,j]
            EstLenAtIntAge[y,s,j] = EstLenAtIntAge[y-1,s,j-1] + EstAnnGrowthInc[y,s,j]
          } # j>1
        } # y>1
      } # j
    } # s
  } # y

  EstLenAtAge = rep(0,nObs)
  tempGrowthInc = rep(0,nObs)
  EstGrowthInc = rep(0,nObs)
  for (j in 1:nObs) {

    y = ObsLenAtAgeDat$YrIndex[j]
    s = ObsLenAtAgeDat$Sex[j]
    i = floor(ObsLenAtAgeDat$Age[j])
    a = ObsLenAtAgeDat$Age[j]

    tempGrowthInc[j] = (Linf[y,s] - EstLenAtIntAge[y,s,i]) * (1 - exp(-vbK[y,s] * (a - i)))
    EstGrowthInc[j] = 1/(1+exp(slope*(bound-tempGrowthInc[j]))) * tempGrowthInc[j]
    EstLenAtAge[j] = EstLenAtIntAge[y,s,i] + EstGrowthInc[j]

  } # j


  SumSq=0
  ObsLenAtAge = ObsLenAtAgeDat$RandLenAtAge
  for (j in 1:nObs) {
    SqRes = (ObsLenAtAge[j] - EstLenAtAge[j]) * (ObsLenAtAge[j] - EstLenAtAge[j])
    SumSq = SumSq + SqRes
  }

  # calculate NLL
  LL = -(nObs / 2.0) * log(SumSq / nObs)
  NLL = -LL
  # cat("NLL",NLL,"Linf_Fem",Linf_Fem[1],"Linf_Mal",Linf_Mal[1],"vbK_Fem",vbK_Fem[1],"vbK_Mal",vbK_Mal[1],"tzero",tzero,'\n')

  EstLinf_Fem <- exp(lnLinf_Fem)
  EstLinf_Mal <- exp(lnLinf_Mal)
  EstvbK_Fem <- exp(lnvbK_Fem)
  EstvbK_Mal <- exp(lnvbK_Mal)
  Esttzero <- tzero

  # RTMB
  RTMB::REPORT(EstLenAtIntAge)
  RTMB::REPORT(ObsLenAtAge)

  RTMB::ADREPORT(EstLinf_Fem)
  RTMB::ADREPORT(EstLinf_Mal)
  RTMB::ADREPORT(EstvbK_Fem)
  RTMB::ADREPORT(EstvbK_Mal)
  RTMB::ADREPORT(Esttzero)

  return(NLL)

}

#' Fit von Betalanffy growth model with year-effects parameters
#'
#' This function fits a von Bertalanffy growth model with year-effects parameters (Linf and/or vbK).
#' The function uses RTMB for optimisation.
#'
#' @param params von Bertalanffy growth parameters (some of which are time-varying)
#'
#' @return SummaryRes (par, convergence, NLL, sdr, params_pt_est, params_sd=params_sd, ParamEst=ParamEst),
#' RTMBReport (EstLenAtIntAge, ObsLenAtAge)
#'
#' @examples
#' library(L3Assess)
#' library(RTMB)
#' set.seed(123)
#' MinAge = 1
#' MaxAge = 40
#' Ages = MinAge:MaxAge
#' NatMort = 0.11
#' FMort = 0.11
#' ZMort = FMort + NatMort
#' SelA50 = c(7,6)
#' SelA95 = c(8,7)
#' SampleSize = 150 # required sample size. For 2 sex model, same sample size generated for each sex.
#' Res=SimAgeFreqData_EqMod(SampleSize, MinAge, MaxAge, SelA50, SelA95, NatMort, FMort)
#' ObsAgeFreq = Res$CatchSample
#' par(mfrow=c(1,1),mar=c(4,3,2,2))
#' plot(MinAge:MaxAge, ObsAgeFreq[1,],col="red",'o')
#' lines(MinAge:MaxAge, ObsAgeFreq[2,],col="blue")
#' # Specify information required to generate length-at-age data with
#' # time-varying growth
#' StartYear = 2000
#' EndYear = 2010
#' nYrs = EndYear - StartYear + 1
#' minIntAge = 0
#' MaxIntAge = 40
#' nSexes = 2
#' # Specify true growth parameter values. Linf and/or vbK can be year-specific.
#' nLinfVals = nYrs
#' nvbKVals = nYrs
#' Linf_Fem = seq(1000,900,-10)
#' Linf_Mal = seq(1100,1000,-10)
#' Linf = data.frame(data.frame(Linf_Fem=Linf_Fem,Linf_Mal=Linf_Mal))
#' vbK_Fem = seq(0.2,0.1,-0.01)
#' vbK_Mal = seq(0.2,0.1,-0.01)
#' vbK = data.frame(data.frame(vbK_Fem=vbK_Fem,vbK_Mal=vbK_Mal))
#' tzero = 0 # Specify common tzero for both sexes to simplify
#' CV = 0.1
#' Res = SimulateTimeVaryingLengthAtAgeData(Linf, vbK, tzero, CV, nLinfVals, nvbKVals, StartYear, EndYear, minIntAge, MaxIntAge, nSexes)
#' ObsLenAtAgeDat = Res$ObsLenAtAgeDat
#' # Fit year-effects growth model in RTMB
#' # Specify starting values
#' InitLinf_Fem = rep(1000,nLinfVals)
#' InitLinf_Mal = rep(1100,nLinfVals)
#' InitvbK_Fem = rep(0.12,nvbKVals)
#' InitvbK_Mal = rep(0.13,nvbKVals)
#' Inittzero = 0
#' # RTMB data input list
#' dat <- list(nLinfVals=nLinfVals,
#'             nvbKVals=nvbKVals,
#'             ObsLenAtAgeDat=ObsLenAtAgeDat)
#' # RTMB param input list
#' params <- list(lnLinf_Fem=log(InitLinf_Fem),
#'                lnLinf_Mal=log(InitLinf_Mal),
#'                lnvbK_Fem=log(InitvbK_Fem),
#'                lnvbK_Mal=log(InitvbK_Mal),
#'                tzero=Inittzero)
#' FittedRes = GetYearEffectGrowthModResults(dat, params)
#' Age = ObsLenAtAgeDat$Age
#' RandLenAtAge = ObsLenAtAgeDat$RandLenAtAge
#' YrIndex = ObsLenAtAgeDat$YrIndex
#' Sex = ObsLenAtAgeDat$Sex
#' IntAge = floor(ObsLenAtAgeDat$Age)
#' # plot simulated length at age data
#' for (y in 1:nYrs) {
#'   for (s in 1:nSexes) {
#'     x=which(YrIndex==y & Sex==s)
#'     if (y==1 & s==1) {
#'       plot(Age[x], RandLenAtAge[x], xlim=c(0,MaxAge), ylim=c(0,1400), col="red")
#'     }
#'     if (y>1 & s==1) {
#'       points(Age[x], RandLenAtAge[x], col="red")
#'     }
#'     if (s==2) {
#'       points(Age[x], RandLenAtAge[x], col="blue")
#'     }
#'   }
#' }
#' # plot female data for some specific ages
#' par(mfcol=c(4,2), mar=c(4,3,2,2))
#' for (s in 1:2) {
#'   for (a in 7:10) {
#'     x=which(Sex==s & IntAge==a)
#'     colour = c("red","blue")
#'     plot(YrIndex[x], RandLenAtAge[x], xlim=c(1,nYrs), ylim=c(300,1200), col=colour[s],main=paste("Age =",a))
#'     lines(1:nYrs,  FittedRes$RTMBReport$EstLenAtIntAge[,s,a])
#'   }
#' }
#' @export
GetYearEffectGrowthModResults <- function(dat, params) {

  # set up for RTMB
  obj <- RTMB::MakeADFun(Calculate_NLL_YearEffectGrowthMod, params)

  # optimising, using nlminb
  nlmb <- nlminb(obj$par, obj$fn, obj$gr)

  # get results
  Res = obj$report()

  RTMBReport <- list(EstLenAtIntAge=Res$EstLenAtIntAge,
                     ObsLenAtAge=Res$ObsLenAtAge)

  # get estimates from ADreport
  sdr <- RTMB::sdreport(obj)
  params_pt_est = as.list(sdr, "Est", report=TRUE)
  params_sd = as.list(sdr, "Std", report=TRUE)

  # store parameter estimates in nicer format
  EstLinf_Fem = data.frame(params_pt_est$EstLinf_Fem, params_pt_est$EstLinf_Fem - 1.96 * params_sd$EstLinf_Fem,
                           params_pt_est$EstLinf_Fem + 1.96 * params_sd$EstLinf_Fem)
  colnames(EstLinf_Fem) = c("Estimate","lw_95%CL","up_95%CL")
  EstLinf_Mal = data.frame(params_pt_est$EstLinf_Mal, params_pt_est$EstLinf_Mal - 1.96 * params_sd$EstLinf_Mal,
                           params_pt_est$EstLinf_Mal + 1.96 * params_sd$EstLinf_Mal)
  colnames(EstLinf_Mal) = c("Estimate","lw_95%CL","up_95%CL")
  EstvbK_Fem = data.frame(params_pt_est$EstvbK_Fem, params_pt_est$EstvbK_Fem - 1.96 * params_sd$EstvbK_Fem,
                          params_pt_est$EstvbK_Fem + 1.96 * params_sd$EstvbK_Fem)
  colnames(EstvbK_Fem) = c("Estimate","lw_95%CL","up_95%CL")
  EstvbK_Mal = data.frame(params_pt_est$EstvbK_Mal, params_pt_est$EstvbK_Mal - 1.96 * params_sd$EstvbK_Mal,
                          params_pt_est$EstvbK_Mal + 1.96 * params_sd$EstvbK_Mal)
  colnames(EstvbK_Mal) = c("Estimate","lw_95%CL","up_95%CL")
  Esttzero = data.frame(params_pt_est$Esttzero, params_pt_est$Esttzero - 1.96 * params_sd$Esttzero,
                        params_pt_est$Esttzero + 1.96 * params_sd$Esttzero)
  colnames(Esttzero) = c("Estimate","lw_95%CL","up_95%CL")

  ParamEst = list(EstLinf_Fem=EstLinf_Fem,
                  EstLinf_Mal=EstLinf_Mal,
                  EstvbK_Fem=EstvbK_Fem,
                  EstvbK_Mal=EstvbK_Mal,
                  Esttzero=Esttzero)


  SummaryRes <- list(nlmb$par,
                     convergence=nlmb$convergence,
                     NLL=nlmb$objective,
                     sdr=sdr,
                     params_pt_est=params_pt_est,
                     params_sd=params_sd,
                     ParamEst=ParamEst)

  results = list(SummaryRes=SummaryRes,
                 RTMBReport=RTMBReport)

  return(results)

}


