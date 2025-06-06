% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{PlotMeanMonthyStGonadWts}
\alias{PlotMeanMonthyStGonadWts}
\title{Plot mean monthly standardised gonad weights}
\usage{
PlotMeanMonthyStGonadWts(
  MatL50,
  SpecFishLength,
  FishRepdat,
  ymax,
  yint,
  GraphTitle,
  xaxis_lab,
  yaxis_lab,
  SampSizelabPosAdj,
  SampSizelab_cex
)
}
\arguments{
\item{MatL50}{length at first maturity}

\item{SpecFishLength}{specified fish lengths to which gonad weights are to be standardised}

\item{FishRepdat}{data frame containing data for individual fish, including month of capture (Mnth), fish length (FishLen),}

\item{ymax}{maximum value for y axis}

\item{yint}{y axis interval}

\item{GraphTitle}{graph title}

\item{xaxis_lab}{x axis label}

\item{yaxis_lab}{y axis label}

\item{SampSizelabPosAdj}{position of monthly sample sizes on plot above upper 95 percent CLs}

\item{SampSizelab_cex}{size of monthly sample size labels}
}
\value{
plot of mean month GSIs with 95 percent confidence limits
}
\description{
This function plots mean monthly standardised gonad weights (for specified fish length value) and associated
95 percent confidence limits, used to assess time and duration of spawning
}
\examples{
# Generate synthetic length composition data for a fished population at equilibrium
library(L3Assess)
set.seed(123)
SampleSize=5000 # sample size for retained catches (and same number for released fish, if an MLL is specified)
MaxAge = 15
TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
NatMort = 4.22/MaxAge
FishMort = NatMort
MaxLen = 500
LenInc = 1
MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA and retention parameters not specified
SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
SelectivityVec = NA # selectivity vector
SelParams = c(150, 20) # L50, L95-L50 for gear selectivity
RetenParams = NA # L50, L95-L50 for retention
DiscMort = 0 # proportion of fish that die due to natural mortality
GrowthCurveType = 1 # 1 = von Bertalanffy, 2 = Schnute
Linf = 350
vbK = 0.4
CVSizeAtAge = 0.08
GrowthParams = c(Linf, vbK)
RefnceAges = NA
Res=SimLenAndAgeFreqData_EqMod(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                         SelParams, RetenParams, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
# plot(Res$ObsDecAgeRetCatch_Fem, Res$ObsRandLenRetCatch_Fem, ylim=c(0,MaxLen), xlim=c(0,MaxAge))
# Now simulate length data from a Dirchlet multinomial distribution, with specified autocorrelation,
# to replicate non-random sampling (i.e. due to limited sampling intensity and fish schooling according to size)
# install.packages("dirmult")
library(dirmult)
set.seed(123)
nSampEvents = 5 # number of sampling events per month
nFishPerSampEvent = 5 # number of fish caught per sampling event
Mnth <- sort(rep(1:12,nSampEvents*nFishPerSampEvent))
nSamples = length(Mnth)
theta_val = 0.9 # level of autocorrelation of lengths within sampling events
midpt = Res$midpt
ExpPropAtLen = Res$ModelDiag$ExpRetCatchPropAtLen_Fem
FishLen <- NA
for (mm in 1:12) {
  res=SimLenFreqDat_DirMultDistn_EqMod(nSampEvents, nFishPerSampEvent, theta_val, midpt, ExpPropAtLen)
  if (mm == 1) {
    FishLen <- rep(res$midpt, res$simLenFreq)
  } else {
    tempFishLen <- rep(res$midpt, res$simLenFreq)
    FishLen <- c(FishLen,tempFishLen)
  }
}
FishLen <- round(FishLen,0)
FishWt <- round(0.00002 * FishLen ^ 3,0)
# set up a linear model with specified coefficients, from which
# fish gonad weight data will be simulated, with error
a = -17.5
b = c(0, -0.2, -0.1, -0.1, -0.1, 0, 0.3, 0.8, 1.1, 1.5, 0.8, 0.2)
c = 3.5
ln_GonadWt <- rep(NA, nSamples)
for (i in 1:nSamples) {
  ln_GonadWt[i] =  (a + b[Mnth[i]] + c*log(FishLen[i])) + rnorm(1,0,1)
}
FishRepdat <- data.frame(GonadWt=exp(ln_GonadWt), ln_GonadWt=ln_GonadWt,
                         FishLen=FishLen, ln_FishLen = log(FishLen),
                         FishWt=FishWt, ln_FishWt = log(FishWt), Mnth)
# inspect simulated data
# par(mfrow = c(3, 2)) # Set the layout of the plots
# plot(FishWt ~ FishLen, data = FishRepdat, main = "Fish weight by Fish length")
# boxplot(FishLen ~ Mnth, data = FishRepdat, main = "FishLen by Month")
# boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_FishLen by Month")
# boxplot(GonadWt ~ Mnth, data = FishRepdat, main = "GonadWt by Month")
# boxplot(ln_FishLen ~ Mnth, data = FishRepdat, main = "ln_GonadWt by Month")
# hist(FishRepdat$FishLen)
# hist(FishRepdat$GonadWt)
# hist(FishRepdat$FishWt)
#
# plot 'true' monthly (relative) trend
par(mfrow=c(2,2))
plot(1:12, b, "o")
# plot GSIs
PlotMeanMonthlyGSIs(MatL50=150, FishRepdat, ymax=30, yint=5, GraphTitle="Females", xaxis_lab="Month",
                    yaxis_lab="GSI", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
# get monthly gonad weights, standardised for length, and plot
res=CalcMeanMonthlyStGonadWts(MatL50=150, SpecFishLength=300, FishRepdat)
res$ModelResults
PlotMeanMonthyStGonadWts(MatL50=150, SpecFishLength=300, FishRepdat, ymax=100, yint=20, GraphTitle="Females",
                         xaxis_lab="Month", yaxis_lab="Stand. gonad wt.", SampSizelabPosAdj=2, SampSizelab_cex=0.8)
}
