\name{spikes}
\alias{spikes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fraud-detection from vote-share data
}
\description{
Implements the resampled kernel density method to detect the excess number of election results with coarse vote-shares (a coarse vote-share is a fraction with a low denominator).
}
\usage{
spikes(data, resamples = 1000, bw = 1e-04, grid = 1001, out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame with three columns with names \eqn{N} (number of votes), \eqn{t} (number who turned out to vote), and \eqn{v} (number who voted for the party the votes of which are being analyzed). Returns error if columns are named incorrectly.
}
  \item{resamples}{
Number of resamples; default \code{1000}
}

  \item{bw}{
Bandwidth for kernel density; default \code{0.0001}
}
  \item{grid}{
Number of points on which the density is estimated; default \code{1001}
}


  \item{out}{
Object containing parameters of beta-mixture model. If \code{spikes} has been called earlier, then \code{out = output$out} will skip density estimation and proceed directly to resampling.
}

}

\value{
\code{spikes} returns object of class \code{out}.
\item{fraud}{Estimated percentage of polling stations with fraud.}
\item{ymax}{Upper envelope of kernel density samples.}
\item{w}{Weights for each bin: the proportion of observations falling into a bin.}
\item{out}{Maximum likelihood estimates of the mixture beta binomial parameters for turnout and votes.}
\item{data}{Data used in estimation.}
}


\seealso{
See Also \code{\link{plot.out}}, \code{\link{summary.out}}
}
\examples{
data(data)
\dontrun{
out <- spikes(data, resamples = 1000)
}
}

