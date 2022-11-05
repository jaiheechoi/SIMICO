\name{sd_off}
\alias{sd_off}
\title{
d^2/d_{theta_k}d{sigma^2}
}
\description{
Calculates the Information matrix term of I_{theta sigma^2} for outcomes l and m.
}
\usage{
sd_off(l, m, phen_l, phen_m, temp_beta, d, apply_diffs, A_i,
   HL_array, HR_array, no_l_all, no_two_all, tpos_all, obs_all, k)
}
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{m}{Index of second outcome of interest.}
  \item{phen_l}{List containing the left and right design matrices and interval times for outcome l.}
  \item{phen_m}{List containing the left and right design matrices and interval times for outcome m.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{d}{Total number of quadrature nodes.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{no_two_all}{n x (k - 2) matrix containing the product of apply_diffs across all outcomes K excluding outcomes l and m.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{k}{Total number of outcomes.}
}
\value{
The output is a matrix containing the component of the information matrix of the sigma and theta parameters.}
