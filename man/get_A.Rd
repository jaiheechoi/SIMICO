\name{get_A}
\alias{get_A}
\title{
Get A vector
}
\description{
Product of difference of survival terms of the left and right interval times, across all outcomes k, summed over all quadrature nodes d.
}
\usage{
get_A(store, weights, d, n)
}
\arguments{
  \item{store}{Matrix of difference of survival values of the left and right time intervals.}
  \item{weights}{Gaussian quadrature weights.}
  \item{d}{Total number of Gaussian quadrature nodes.}
  \item{n}{Total number of observations.}
}
\value{
The output is a vector used to compute the derivative terms.
}
