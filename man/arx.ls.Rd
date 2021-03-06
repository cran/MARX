% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MARX_functions.R
\name{arx.ls}
\alias{arx.ls}
\title{The ARX estimation by OLS function}
\usage{
arx.ls(y, x, p)
}
\arguments{
\item{y}{Data vector of time series observations.}

\item{x}{Matrix of data (every column represents one time series). Specify NULL or "not" if not wanted.}

\item{p}{Number of autoregressive terms to be included.}
}
\value{
\item{coefficients}{Vector of estimated coefficients.}

\item{coef.auto}{Vector of estimated autoregressive parameters.}

\item{coef.exo}{Vector of estimated exogenous parameters.}

\item{mse}{Mean squared error.}

\item{residuals}{Residuals.}

\item{loglikelihood}{Value of the loglikelihood.}

\item{fitted.values}{Fitted values.}

\item{df}{Degrees of freedom.}

\item{vcov}{Variance-covariance matrix of residuals.}
}
\description{
This function allows you to estimate ARX models by ordinary least squares (OLS).
}
\examples{
data <- sim.marx(c('t',3,1),c('t',1,1),100,0.5,0.4,0.3)
arx.ls(data$y,data$x,2)
}
\author{
Sean Telg
}
\keyword{estimation}
\keyword{pseudo-causal}
