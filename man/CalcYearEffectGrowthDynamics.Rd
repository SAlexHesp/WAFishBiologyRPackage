% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{CalcYearEffectGrowthDynamics}
\alias{CalcYearEffectGrowthDynamics}
\title{Simulate multiple years of length at age data with time-varying growth}
\usage{
CalcYearEffectGrowthDynamics(ObsLenAtAgeDat, Linf, vbK, tzero)
}
\arguments{
\item{ObsLenAtAgeDat}{year specific length-at-age data}

\item{Linf}{asymptotic length (by sex, and year, if multiple years)}

\item{vbK}{growth coefficient (by sex, and year, if multiple years)}

\item{tzero}{combined sex, single value}
}
\value{
results (EstLenAtIntAge, EstLenAtAge)
}
\description{
This function simulates multiple connsecutive years of length-at-age data
with time-varying growth, based on year specific von Bertalanffy growth
parameters (Linf and/or vbK). Data are simulated assuming constant mortality
and selectivity.
}
\keyword{internal}
