% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_Wilmoth.R
\name{find_my_case}
\alias{find_my_case}
\title{Function that determines the case/problem we have to solve
It also performs some checks}
\usage{
find_my_case(q0_5, q0_1, q15_45, q15_35, e0, k)
}
\arguments{
\item{q0_5}{5q0. The probability that a new-born will die during the 
subsequent 5 years;}

\item{q0_1}{1q0. The probability that a life aged 0 will die during the 
following year;}

\item{q15_45}{45q15. The probability that a life aged 15 will die during 
the subsequent 45 years;}

\item{q15_35}{35q15. The probability that a life aged 15 will die during 
the subsequent 35 years;}

\item{e0}{Life expectancy at birth;}

\item{k}{k-parameter in the log-quadratic model;}
}
\description{
Function that determines the case/problem we have to solve
It also performs some checks
}
\keyword{internal}
