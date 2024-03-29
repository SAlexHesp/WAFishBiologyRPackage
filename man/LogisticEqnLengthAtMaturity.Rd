% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{LogisticEqnLengthAtMaturity}
\alias{LogisticEqnLengthAtMaturity}
\title{Logistic curve for length at maturity}
\usage{
LogisticEqnLengthAtMaturity(params)
}
\arguments{
\item{params}{c(L50, L95) or c(Pmax, L50, L95) for CurveType=1, or c(Q, B, V) or c(Pmax, Q, B, V) for CurveType=2}
}
\value{
probability of maturity for observed lengths. Used in calculation of objective function
for length at maturity model
}
\description{
This function applies a logistic curve (symmetric or asymmetric) for describing probability of maturity at length
}
\keyword{internal}
