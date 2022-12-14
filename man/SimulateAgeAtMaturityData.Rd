% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{SimulateAgeAtMaturityData}
\alias{SimulateAgeAtMaturityData}
\title{Simulate age at maturity data}
\usage{
SimulateAgeAtMaturityData(nSamples, nSexes, MinAge, MaxAge, MaturityParams)
}
\arguments{
\item{nSamples}{number of required samples}

\item{nSexes}{number of sexes}

\item{MinAge}{minimum age}

\item{MaxAge}{maximum age}

\item{A50}{age at 50 percent maturity}

\item{A95}{age at 95 percent maturity}
}
\description{
Simulate age at maturity data, given specified growth, size at maturity parameters and required sample size
}
\examples{
# generate some synthetic maturity data
set.seed(123)
MinAge = 0
MaxAge = 20
nSexes = 1
nSamples = 300
A50 = 4
A95 = 6
Pmax = 1.0
MaturityParams = c(A50, A95, Pmax)
# nSexes = 2
# nSamples = c(300,300)
# A50 = c(4,4.5)
# A95 = c(6,6)
# Pmax = c(1.0,1.0)
# MaturityParams = c(A50, A95, Pmax)
SimulateAgeAtMaturityData(nSamples, nSexes, MinAge, MaxAge, MaturityParams)
}
