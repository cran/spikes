\name{summary.out}
\alias{summary.out}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize
}
\description{
Extracts estimate of fraud and 95 percent credible interval (if such is estimated) for the object of class \code{out} returned by \code{spikes} or \code{confint.out}.
}
\usage{
\method{summary}{out}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{out}}
  \item{...}{additional arguments}
}


\note{
If the argument is from \code{spikes}, then \code{summary} returns degenerate credible interval, as it was not estimated. To return a proper credible interval, \code{confint.out} must take as its argument object returned by \code{confint.out}.
}

\examples{
data(output)
summary(output)
}
