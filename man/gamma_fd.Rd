\name{gamma_fd}
\alias{gamma_fd}
\title{
d/d_{gamma_l}
}
\description{
Calculates the gradient term for U_g for the score statistic.
}
\usage{
gamma_fd(l, HL_array, HR_array, tpos_all, obs_all,
   temp_beta, A_i, no_l_all, gMat, a1, a2, d)
}
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{temp_beta}{Vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{gMat}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{d}{Number of quadrature nodes.}
}
\value{
The output is a vector containing the first derivative with respect to gamma.
}
