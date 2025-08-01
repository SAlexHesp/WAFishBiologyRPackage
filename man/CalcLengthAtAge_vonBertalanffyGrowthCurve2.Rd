% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{CalcLengthAtAge_vonBertalanffyGrowthCurve2}
\alias{CalcLengthAtAge_vonBertalanffyGrowthCurve2}
\title{Calculate expected lengths at age from von Bertalanffy growth curve}
\usage{
CalcLengthAtAge_vonBertalanffyGrowthCurve2(
  params,
  DataType,
  nSexes,
  ObsAge,
  plotages
)
}
\arguments{
\item{params}{c(log(Linf),log(vbK),tzero)}

\item{DataType}{# 1=lengths at age data for individual fish (single sex, or sexes recorded),
2=mean length at age and sd data from mixture analysis, 3=lengths at age data for individual fish
(two sex, sexes not recorded)}

\item{nSexes}{number of sexes}

\item{ObsAge}{observed ages}

\item{plotages}{specified ages for plotting}
}
\value{
specified ages (plotages) and expected lengths at ages (plotlengths)
}
\description{
Calculates expected lengths at age from von Bertalanffy growth curve
and associated growth curve parameter values. Used in plotting where
ages are pres-specified
}
\keyword{internal}
