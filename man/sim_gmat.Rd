\name{sim_gmat}
\alias{sim_gmat}
\title{
Simulate genetic matrix.
}
\description{
Simulates a n x q genetic matrix with the option to specify the common pairwise correlation.
}
\usage{
sim_gmat(n,q,rho)
}
\arguments{
  \item{n}{Total number of observations.}
  \item{q}{Total number of SNPs.}
  \item{rho}{Common pairwise correlation parameter.}
}
\value{
The result of a n x q genetic matrix of q SNPs.
}
\examples{
# Set sample size
n = 100

# Set number of SNPs
q = 5

# Set common pairwise correlation
rho = 0.1

# Simulate genetic matrix
gMat <- sim_gmat(n, q, rho)
}
