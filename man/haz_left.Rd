\name{haz_left}
\alias{haz_left}
\title{
H_{ik}(L_{ik})
}
\description{
Calculates the hazard function of the left time interval for outcome l.
}
\usage{
haz_left(l, d, temp_beta, phen, r1, k)
}
\arguments{
  \item{l}{Outcome of interest.}
  \item{d}{Total number of Gaussian quadrature nodes.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{phen}{list of data matrices containing both left and right information.}
  \item{r1}{Gaussian quadrature nodes.}
  \item{k}{Total number of outcomes.}
}
\value{
The output is a vector of the hazard values of the left times.
}
