\name{ss_fd}
\alias{ss_fd}
\title{
d/d_{sigma^2}
}
\description{
Calculates the first derivative term with respect to sigma^2.}
\usage{
ss_fd(l, phen, HL_array, HR_array, tpos_all, obs_all,
   apply_diffs, temp_beta, A_i, no_l_all, k, d)
}
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{phen}{List containing all the left and right design matrices.}
  \item{HL_array}{n x k matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x k matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes k summed over all quadrature nodes d.}
  \item{no_l_all}{n x (k - 1) matrix containing the product of apply_diffs across all outcomes k excluding the current outcome l.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}
}
\value{
The output is a single value for the first derivative with respect to sigma.}
