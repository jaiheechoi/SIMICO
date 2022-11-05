\name{gamma_off}
\alias{gamma_off}
\title{
d^2/d_{gamma_l}d{gamma_m}
}
\description{
Calculates the [off-diagonal] Information matrix term for I_{gamma gamma} with respect to outcome l and outcome m.
}
\usage{
gamma_off(l, m, HL_array, HR_array,
   tpos_all, obs_all, temp_beta, A_i,
   no_l_all, no_two_all, gMat, a1, a2, k, d)
}
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{m}{Index of second outcome of interest.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{temp_beta}{Vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{no_two_all}{n x (k - 2) matrix containing the product of apply_diffs across all outcomes k excluding the outcomes l and m.}
  \item{gMat}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}
}
\value{
The output is a matrix containing the component of the information matrix of the gamma parameter for outcomes l and m.
}
