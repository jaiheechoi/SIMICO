\name{st_off}
\alias{st_off}
\title{
d^2/d_{theta_l}d{sigma^2}
}
\description{
Calculates the Information matrix term of I_{eta theta} for one outcome of interest l.
}
\usage{
st_off(l, HL_array, HR_array, xAll, apply_diffs,
   temp_beta, A_i, no_l_all, no_two_all, k, d)
}
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{xAll}{List containing the left and right matrices and event times.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{no_two_all}{n x (k - 2) matrix containing the product of apply_diffs across all outcomes k excluding outcomes l and m.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}
}
\value{
The output is a matrix containing the component of the information matrix of the theta eta parameters for outcome l.}
