\name{plot.out}
\alias{plot.out}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots output of \code{spikes}
}
\description{
Plots the observed kernel density of data and the upper envelope of the resampled densities.
}
\usage{
\method{plot}{out}(x, main = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class \code{out}
}
  \item{main}{
Title, \code{NULL} by default
}
  \item{...}{
additional plotting arguments
}
}
\examples{
data(output)
plot(output)
}
