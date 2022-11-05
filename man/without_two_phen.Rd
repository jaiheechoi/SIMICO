\name{without_two_phen}
\alias{without_two_phen}
\title{
Survival Difference Product without Outcomes l and m
}
\description{
Differnence of survival functions multiplied across all outcomes excluding outcomes l and m.
}
\usage{
without_two_phen(l, m, k, store, n, d)
}
\arguments{
  \item{l}{The first outcome of interest.}
  \item{m}{The second outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{store}{Array of difference between left and right survival values.}
  \item{n}{Total number of observation.}
  \item{d}{Total number of quadrature nodes.}
}
\value{
A n x (k-2) matrix containing the product of all the differences of left and right survival values across all outcomes excluding outcomes l and m.
}
