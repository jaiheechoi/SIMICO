\name{fd_term}
\alias{fd_term}
\title{
d/d_{theta_l}
}
\description{
Calculate the first derivative of the theta terms for outcome l.
}
\usage{
fd_term(l, temp_beta, phen,d, apply_diffs,
   A_i, no_l_all,HL_array, HR_array)
}
\arguments{
  \item{l}{Outcome of interest.}
  \item{temp_beta}{Vector of fitted coefficients.}
  \item{phen}{list containing the covariate design matrices.}
  \item{d}{Number of quadrature nodes.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
}
\value{
The output is a 1 x (p + 2) vector of the first derivative terms for outcome l.
}
