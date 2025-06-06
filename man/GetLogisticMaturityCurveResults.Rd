% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{GetLogisticMaturityCurveResults}
\alias{GetLogisticMaturityCurveResults}
\title{Get outputs from a fitted logistic maturity curve.}
\usage{
GetLogisticMaturityCurveResults(
  params,
  nSexes,
  LogisticModType,
  CurveType,
  ObsLen,
  ObsAgeCl,
  ObsMatCat,
  ErrOpt
)
}
\arguments{
\item{params}{c(L50,L95) length at 50 and 95 percent maturity (L50)}

\item{nSexes}{c(L50,L95) number of sexes}

\item{LogisticModType}{1=length-based, 2=age-based}

\item{CurveType}{1=symmetric logistic, 2=asymmetric logistic}

\item{ObsLen}{vector of observed lengths (set to NA for age model)}

\item{ObsAgeCl}{vector of observed age classes (set to NA for length model)}

\item{ObsMatCat}{vector of observed maturity categories (0=immature, 1=mature)}

\item{ErrOpt}{method for uncertainty calculation (1=varcov approx, 2=bootstrap)}
}
\value{
negative log-likelihood (nll), nlminb convergence diagnostic (convergence), sample size (SampleSize)
maturity parameter estimates with lower and upper 95 percent confidence limits (ParamEst), point estimates
for growth parameters (params), variance-covariance matrix (vcov.params)
}
\description{
This function fits a logistic maturity curve to a sample of maturity-at-length data
by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
It provides various statistical outputs in include convergence statistics,
parameter estimated and associated 95 percent confidence limits and associated variance-covariance matrix,
calculated using the MASS package.
}
\examples{
# generate synthetic length at maturity data
set.seed(123)
MaxAge = 20
MinLen = 0
MaxLen = 400
LenInc = 20
# single sex/combined sexes
nSexes = 1
nSamples = 300
CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
Linf = 300
vbK = 0.3
tzero = 0
CVSizeAtAge = 0.1
GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
L50 = 220
L95 = 270
Pmax = 1.0
MaturityParams = c(L50, L95, Pmax)
res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
# fit model to synthetic data
ObsAgeCl=NA
ObsLen=res$ObsLen
ObsMatCat=res$ObsMatCat
LogisticModType = 1 # 1=length, 2=age
CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
ErrOpt = 1 # 1=varcov approx, 2=bootstrap
# 2 parameter model symmetric curve
InitL50 = 200
InitL95 = 250
params = c(InitL50, InitL95) # without Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
# 3 parameter model symmetric curve
InitL50 = 200
InitL95 = 250
InitPmax = 0.9
InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
params = c(InitPmax_logit, InitL50, InitL95) # with Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
# 3 parameter model asymmetric curve
CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
Q = 20
B = 200
V = 2
Pmax = 1.0
MaturityParams = c(Q, B, V, Pmax)
res = SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
plot(res$midpt, res$PropMat)
# fit model to synthetic data
ObsAgeCl=NA
ObsLen=res$ObsLen
ObsMatCat=res$ObsMatCat
LogisticModType = 1 # 1=length, 2=age
ErrOpt = 1 # 1=varcov approx, 2=bootstrap
InitQ = 25
InitB = 250
InitV = 1.5
params = c(InitQ, InitB, InitV) # with Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
# separate sexes
nSexes = 2
nSamples = c(500,500)
CurveType = 1 # 1 = symmetric logistic, 2 = asymmetric logistic
Linf = c(300,350)
vbK = c(0.3,0.3)
tzero = c(0,0)
CVSizeAtAge = c(0.1,0.1)
GrowthParams = c(Linf, vbK, tzero, CVSizeAtAge)
L50 = c(200,220)
L95 = c(240,250)
Pmax = c(1,1)
MaturityParams = c(L50, L95, Pmax)
res=SimulateLengthAtMaturityData(nSamples, CurveType, nSexes, MaxAge, MinLen, MaxLen, LenInc, GrowthParams, MaturityParams)
ObsAgeCl=NA
FemObsLen=res$FemObsLen
MalObsLen=res$MalObsLen
ObsLen = as.matrix(t(data.frame(FemObsLen=FemObsLen,MalObsLen=MalObsLen)))
FemObsMatCat=res$FemObsMatCat
MalObsMatCat=res$MalObsMatCat
ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
LogisticModType = 1 # 1=length, 2=age
ErrOpt = 1 # 1=varcov approx, 2=bootstrap
# 2 parameter model symmetric curve
InitL50 = c(200, 220)
InitL95 = c(250, 270)
params = c(InitL50, InitL95) # without Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
# 3 parameter model symmetric curve
InitL50 = c(200, 220)
InitL95 = c(250, 270)
InitPmax = c(0.8, 0.8)
InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
params = c(InitPmax_logit, InitL50, InitL95) # with Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get length-at-maturity results
# age-at-maturity
MinAge = 0
MaxAge = 20
nSexes = 1
CurveType = 1
nSamples = 300
A50 = 4
A95 = 6
Pmax = 1.0
MaturityParams = c(A50, A95, Pmax)
res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
ObsLen=NA
ObsAgeCl=res$ObsAgeCl
ObsMatCat=res$ObsMatCat
LogisticModType = 2 # 1=length, 2=age
ErrOpt = 1 # 1=varcov approx, 2=bootstrap
# 2 parameter model
InitA50 = 3
InitA95 = 5
params = c(InitA50, InitA95) # without Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
# 3 parameter model
InitA50 = 4.5
InitA95 = 7
InitPmax = 0.9
InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
params = c(InitPmax_logit, InitA50, InitA95) # with Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
# # 3 parameter age-based asymmetric curve
# # generate synthetic length at maturity data
set.seed(144)
MinAge = 0
MaxAge = 20
nSexes = 1
nSamples = 500
CurveType = 2 # 1 = symmetric logistic, 2 = asymmetric logistic
Q = 1 # controls spread
B = 5 # inflection point
V = 2 # controls skew or asymmetry
Pmax = 1
MaturityParams = c(Q,B,V,Pmax)
res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
ObsLen=NA
ObsAgeCl=res$ObsAgeCl
ObsMatCat=res$ObsMatCat
LogisticModType = 2 # 1=length, 2=age
ErrOpt = 1 # 1=varcov approx, 2=bootstrap
InitQ = 2 # related to y(0) value
InitB = 6 # the growth rate
InitV = 3 # affects near which asymptote maximum growth occurs
params = c(log(InitQ), InitB, log(InitV)) # without Pmax parameter
Res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
# separate sexes
nSexes = 2
nSamples = c(300,300)
CurveType = 1
A50 = c(4,4.5)
A95 = c(6,6)
Pmax = c(0.95,0.95)
MaturityParams = c(A50, A95, Pmax)
res=SimulateAgeAtMaturityData(nSamples, CurveType, nSexes, MinAge, MaxAge, MaturityParams)
FemObsAgeCl=res$FemObsAgeCl
MalObsAgeCl=res$MalObsAgeCl
ObsAgeCl = as.matrix(t(data.frame(FemObsAgeCl=FemObsAgeCl,MalObsAgeCl=MalObsAgeCl)))
FemObsMatCat=res$FemObsMatCat
MalObsMatCat=res$MalObsMatCat
ObsMatCat = as.matrix(t(data.frame(FemObsMatCat=FemObsMatCat,MalObsMatCat=MalObsMatCat)))
# 2 parameter model
InitA50 = c(3,3)
InitA95 = c(5,5)
params = c(InitA50, InitA95) # without Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
# 3 parameter model
InitA50 = c(4.5,4)
InitA95 = c(7,7)
InitPmax = c(0.9,0.9)
InitPmax_logit = log(InitPmax/(1-InitPmax)) # logit transform
params = c(InitPmax_logit, InitA50, InitA95) # with Pmax parameter
res=GetLogisticMaturityCurveResults(params, nSexes, LogisticModType, CurveType, ObsLen, ObsAgeCl, ObsMatCat, ErrOpt) # get age-at-maturity results
}
