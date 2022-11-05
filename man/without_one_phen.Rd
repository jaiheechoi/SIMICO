\name{without_one_phen}
\alias{without_one_phen}
\title{
Survival Difference Product without Outcome l
}
\description{
Calculate the product of the difference between survival terms excluding that of the outcome of interest.
}
\usage{
without_one_phen(l, k, store)
}
\arguments{
  \item{l}{Outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{store}{Array of difference between left and right survival values.}
}
\value{
A n x (k-1) matrix where each column is the product of all the differences of left and right survival values across all outcomes excluding the column index outcome.
}
