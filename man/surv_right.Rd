\name{surv_right}
\alias{surv_right}
\title{
S_{ik}(R_{ik})
}
\description{
Calculates the survival function of the right time interval for outcome k for quadrature node d.
}
\usage{
surv_right(l, d, temp_beta, phen, r1, k)
}
\arguments{
  \item{l}{Outcome of interest.}
  \item{d}{Total number of Gaussian quadrature nodes.}
  \item{temp_beta}{Vector of fitted coefficients.}
  \item{phen}{List of data matrices containing both left and right information.}
  \item{r1}{Gaussian quadrature nodes.}
  \item{k}{Total number of outcomes.}
}
\value{
The output is a vector of the survival values of the left times.
}
