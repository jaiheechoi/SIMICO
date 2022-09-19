\name{simico_fit_null}
\alias{simico_fit_null}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simico_fit_null()
}
\description{
Fit the null model via newton raphson for multiple outcomes interval-censored skat.
}
\usage{
simico_fit_null(init_beta, epsilon, xDats, lt_all, rt_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta_fit}{Starting values for NR.}
  \item{epsilon}{Stopping criterion for NR.}
  \item{xDats}{List of left and right design matrices.}
  \item{lt_all}{n x k matrix of left times.}
  \item{rt_all}{n x k matrix of right times.}
  \item{k}{Total number of outcomes.}
  \item{d}{Total number of quadrature nodes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  temp_beta = temp_beta, iter = iter, diff = diff, jmat = jmat, grad = grad
 \item{beta_fit}{vector of fitted coefficients}
 \item{iter}{number of iternations needed for the Newton Raphson to converge}
 \item{diff}{Difference between the current values of temp_beta and the previous iteration of temp_beta}
 \item{jmat}{Information matrix of the theta parameters}
 \item{grad}{Vector of the first derivaive of the theta parameters}
%% ...
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

# Set number of quadrature nodes
d = 100

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

# Set the initial estimate values
init_beta <-c (rep(c(0, 0, 0, 1, 0), k), 1)

# Run the newton-raphson
nullFit <- simico_fit_null(init_beta = init_beta, epsilon = 10^-5, xDats = exampleDat$sample$xDats, lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, k = k, d = d)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line