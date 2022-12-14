% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WA_biology_methods.R
\name{MixtureDistnCurve_SizeFreq}
\alias{MixtureDistnCurve_SizeFreq}
\title{Get expected size frequency from size mixture distribution model}
\usage{
MixtureDistnCurve_SizeFreq(params)
}
\arguments{
\item{params}{c(Mean1, sd1) normal distn: 1 mode, or
c(Shape1, Rate1) gamma distn: 1 mode, or
c(Mean1, Mean2, sd1, PropZero) normal distn: 2, common sd, or
c(Shape1, Shape2, Rate1, PropZero) gamma distn: 2, common sd, or
c(Mean1, Mean2, sd1, sd2, PropZero) normal distn: 2, separate sds, or
c(Shape1, Shape2, Rate1, Rate2, PropZero) gamma distn: 2 modes, separate sds, or
c(Mean1, Mean2, Mean3, sd1, sd2, sd3, PropZero, PropOne) normal distn: 3 modes, separate sds, or
c(Shape1, Shape2, Shape3, Rate1, Rate2, Rate3, PropZero, PropOne) gamma distn: 3 modes, separate sds}
}
\value{
Expected size frequency at length (ExpFreqAtLen)
}
\description{
Gets expected size frequency from size mixture distribution model, used in calculating
confidence limits for mixture distribution curves applying a parametric resampling procedure
}
\keyword{internal}
