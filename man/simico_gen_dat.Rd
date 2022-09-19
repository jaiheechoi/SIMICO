\name{simico_gen_dat}
\alias{simico_gen_dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simico_gen_dat()
}
\description{
Generate multiple interval-censored data under proportional hazards model.

}
\usage{
simico_gen_dat(bhFunInv, obsTimes, windowHalf = 0.1, n, k, tauSq, gMatCausal, effectSizes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
% bhFunInv, obsTimes = 1:3, windowHalf = 0.1, n, p, k, tauSq, gMatCausal, effectSizes
  \item {bhFunInv}{The inverse of the baseline hazard function.}
  \item{obsTimes}{Vector of the intended observation times.}
  \item{windowHalf}{The amount of time before or after the intended obsTimes that a visit might take place.}
  \item{n}{Total number of observations.}
  \item{p}{Total number of covariates.}
  \item{k}{Total number of outcomes.}
  \item{tausq}{Variance of the subject specific random effect.}
  \item{gMatCausal}{Matrix of subsetted genetic information for only a select causal SNPs.}
  \item{effectSizes}{Vector of genetic effect sizes. Should be entered as a vector the same length as the number of outcomes. Each effect size will be multiplied by the }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
%rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd
%%  If it is a LIST, use
% exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
              rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd, dataFull = xAll))

%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
   \item{exactTimesMat}{n x k matrix containing the simulated exact times that the event occurred}
   \item{leftTimesMat}{n x k matrix containing the left time interval that is observed}
   \item{rightTimesMat}{n x k matrix containing the right time interval that is observed}
   \item{obsInd}{n x k matrix containing a indictor for whether that time is right-censored or not}
   \item{tposInd}{n x k matrix containing a indictor for whether that time is left-censored or not}
   \item{fullDat}{Data in complete form to enter into SIMICO functions.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Set three outcomes
k = 3

# Set number of observations
n = 1500

# Set number of covariates
p = 2

# Set number of SNPs
q = 50

# Variance of subject-specific random effect
tausq = 1

# Define the effect sizes
effectSizes <- c(.03, .15, .3)

# the baseline cumulative hazard function
bhFunInv <- function(x) {x}

set.seed(1)
# Fixed effects
xMat <- cbind(rnorm(n), rbinom(n=n, size=1, prob=0.5))

# Genetic effects
gMat <- sim_gmat(n, q, rho)

# True model has nothing
fixedMat <- matrix(data=0, nrow=n, ncol=k)

# Generate the multiple outcomes
exampleDat <- newgen(bhFunInv = bhFunInv, obsTimes = 1:3, windowHalf = 0.1, n, p, k, tauSq, gMat, effectSizes)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line