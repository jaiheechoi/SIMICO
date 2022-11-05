\name{simico_gen_dat}
\alias{simico_gen_dat}
\title{
simico_gen_dat()
}
\description{
Generate multiple interval-censored data under proportional hazards model.

}
\usage{
simico_gen_dat(bhFunInv, obsTimes = 1:3, windowHalf = 0.1,
   n, p, k, tauSq, gMatCausal, xMat, effectSizes)
}
\arguments{
  \item{bhFunInv}{The inverse of the baseline hazard function.}
  \item{obsTimes}{Vector of the intended observation times.}
  \item{windowHalf}{The amount of time before or after the intended obsTimes that a visit might take place.}
  \item{n}{Total number of observations.}
  \item{p}{Total number of covariates.}
  \item{k}{Total number of outcomes.}
  \item{tauSq}{Variance of the subject specific random effect.}
  \item{gMatCausal}{Matrix of subsetted genetic information for only a select causal SNPs.}
  \item{xMat}{Matrix of covariates.}
  \item{effectSizes}{Vector of genetic effect sizes. Should be entered as a vector the same length as the number of outcomes.}
}
\value{
   \item{exactTimesMat}{n x k matrix containing the simulated exact times that the event occurred.}
   \item{leftTimesMat}{n x k matrix containing the left time interval that is observed.}
   \item{rightTimesMat}{n x k matrix containing the right time interval that is observed.}
   \item{obsInd}{n x k matrix containing a indictor for whether that time is right-censored or not.}
   \item{tposInd}{n x k matrix containing a indictor for whether that time is left-censored or not.}
   \item{fullDat}{Data in complete form to enter into SIMICO functions.}
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
}
