\name{confInt}
\alias{confInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Credible interval
}
\description{
Estimates credible interval by (1) taking a draw from the posterior density (2) implementing the RKD step. The procedure is repeated \code{boot} times.
}
\usage{
confInt(object, boots = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{object of class \code{out} returned by \code{spikes}}
  \item{boots}{number of samples from the posterior; defaul \code{100}}
}
