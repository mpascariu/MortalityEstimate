% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_Wilmoth.R
\name{wilmoth}
\alias{wilmoth}
\title{Fit the Log-Quadratic Model}
\usage{
wilmoth(x, mx = NULL, LT = NULL, verbose = TRUE, control = list())
}
\arguments{
\item{x}{Numerical vector containing ages corresponding to the input data 
(in 'mx' or 'LT').}

\item{mx}{Input data. A data.frame / matrix with death rates.}

\item{LT}{Input data. A collection of life tables. If \code{mx} is provided 
\code{LT} is not necessary and vice versa.}

\item{verbose}{Logical. Choose whether to display a progress bar during the 
fitting process. Default: TRUE.}

\item{control}{List with additional parameters. 
See \code{\link{wilmoth_control}}.}
}
\value{
The output is of class \code{wilmoth} with the components:
 \item{input}{A list with input objects provided in the function;}
 \item{call}{An unevaluated function call, that is, an unevaluated 
 expression which consists of the named function applied to the given 
 arguments;}
 \item{coefficients}{Estimated coefficients;}
 \item{k}{The estimated k parameter;}
 \item{fitted.values}{Fitted values of the estimated model;}
 \item{residuals}{Deviance residuals;}
 \item{model.info}{ Model details (equation). The relationship between 
the death rate at age x, and the probability of dying between birth and 
age 5.}
}
\description{
Estimate the Log-Quadratic model. The implemented estimation using the 
bi-weight procedure is described in the Appendix of 
\insertCite{wilmoth2012;textual}{MortalityEstimate}
}
\examples{
\dontrun{
# DATA
HMD719f <- HMD719[HMD719$sex == "female", ]

# Fit model
x <- c(0,1, seq(5, 110, by = 5))
W <- wilmoth(x, LT = HMD719f)
}
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{wilmothLT}}
}
