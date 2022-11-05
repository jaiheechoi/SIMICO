\name{simico_fit_null}
\alias{simico_fit_null}
\title{
simico_fit_null()
}
\description{
Fit the null model via newton raphson for multiple outcomes interval-censored skat.
}
\usage{
simico_fit_null(init_beta, epsilon, xDats, lt_all, rt_all, k, d)
}
\arguments{
  \item{init_beta}{Starting values for NR.}
  \item{epsilon}{Stopping criterion for NR.}
  \item{xDats}{List of left and right design matrices.}
  \item{lt_all}{n x k matrix of left times.}
  \item{rt_all}{n x k matrix of right times.}
  \item{k}{Total number of outcomes.}
  \item{d}{Total number of quadrature nodes.}
}
\value{
 \item{beta_fit}{Vector of fitted coefficients.}
 \item{iter}{Number of iterations needed for the Newton-Raphson to converge.}
 \item{diff}{Difference between the current values of temp_beta and the previous iteration of temp_beta.}
 \item{jmat}{Information matrix of the theta parameters.}
 \item{grad}{Vector of the first derivaive of the theta parameters.}
}
\examples{
# Set number of outcomes
k = 2

# Set number of observations
n = 100

# Set number of covariates
p = 2

# Set number of SNPs
q = 50

# Set number of causal SNPs
num = 5

# Set number of quadrature nodes
d = 100

# Variance of subject-specific random effect
tauSq = 1

# Define the effect sizes
effectSizes <- c(.03, .15)

# Set MAF cutoff for causal SNPs
Causal.MAF.Cutoff = 0.1

# the baseline cumulative hazard function
bhFunInv <- function(x) {x}

set.seed(1)

# Generate covariate matrix
xMat <- cbind(rnorm(n), rbinom(n=n, size=1, prob=0.5))

# Generate genetic matrix
gMat <- matrix(data=rbinom(n=n*q, size=2, prob=0.1), nrow=n)

# Get indices to specific select causal variants
idx <- Get_CausalSNPs_bynum(gMat, num, Causal.MAF.Cutoff)

# Subset the gMat
gMatCausal <- gMat[,idx]

# Generate the multiple outcomes
exampleDat <- simico_gen_dat(bhFunInv = bhFunInv, obsTimes = 1:3,
                             windowHalf = 0.1, n, p, k, tauSq, gMatCausal,
                             xMat, effectSizes)

# Set the initial estimate values
init_beta <-c (rep(c(0, 0, 0, 1, 0), k), 1)

# Run the Newton-Raphson
nullFit <- simico_fit_null(init_beta = init_beta,
   epsilon = 10^-5, xDats = exampleDat$fullDat$xDats,
   lt_all = exampleDat$leftTimesMat,
   rt_all = exampleDat$rightTimesMat,
   k = k, d = d)

}
