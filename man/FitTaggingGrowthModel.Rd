% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{FitTaggingGrowthModel}
\alias{FitTaggingGrowthModel}
\title{Fit a tagging growth model to tag-recapture data}
\usage{
FitTaggingGrowthModel(
  params,
  nstep,
  Obs_delta_t,
  Obs_Initlen,
  Obs_Finlen,
  nobs
)
}
\arguments{
\item{params}{log(c(L50_1, L95_1, L50_2, L95_2, Max_increment)) double logistic model, or
log(c(Gaussian_A, Gaussian_u, Gaussian_sd)) Gaussian function, or log(c(vb_Linf, vb_K)) von Bertalanffy, or
log(c(Gomp_Linf, Gomp_G)) Gompertz}

\item{nstep}{number of numerical integration steps (higher number increases accuracy but reduces program speed)}

\item{Obs_delta_t}{observed durations at liberty for individual animals}

\item{Obs_Initlen}{observed initial lengths for individual animals}

\item{Obs_Finlen}{observed final lengths for individual animals}

\item{nobs}{number of observations}
}
\value{
nlmb (stored output from internal R nlminb optimisation function)
}
\description{
This function fits a tagging growth model to a sample of fish tag-recapture data
by minimising the negative log-likelihood associated with the parameters and data, using nlminb.
}
\keyword{(Numeric}
\keyword{Obs_Finlen}
\keyword{Obs_Initlen,}
\keyword{Obs_delta_t,}
\keyword{Vectors)}
\keyword{internal}
