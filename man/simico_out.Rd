\name{simico_out}
\alias{simico_out}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get P-Values
}
\description{
Calculate test statistic and p-values for multiple outcome test and multiple burden test.
}
\usage{
simico_out(nullFit, xDats, lt_all, rt_all, Itt, a1, a2, G, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
% nullFit, xDats, lt_all, rt_all, Itt, a1, a2, G, k, d
  \item{nullFit}{Results of the Newton-Raphson: estimates of the beta coefficients.}
  \item{xDats}{List of design matrices.}
  \item{lt_all}{Matrix containing the generated left interval times.}
  \item{rt_all}{Matrix containing the generated right interval times.}
  \item{Itt}{I_{\theta\theta} - Information matrix of thetha.}
  \item{G}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  multQ = gamma_score, multP = pval$Qq, burdQ = burden_score, burdP = p_burden
  \item{multQ}{Score statistic for multiple outcome test.}
  \item{multP}{P-value for multiple outcome test.}
  \item{burdQ}{Score statistic for multiple burden test.}
  \item{burdP}{P-value for multiple burden test.}
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

# Get indices to specific select causal variants
idx <- Get_CausalSNPs_bynum(gMat, num, Causal.MAF.Cutoff)

# Subset the gMat
gMatCausal <- gMat[,idx]

# True model has nothing
fixedMat <- matrix(data=0, nrow=n, ncol=k)

# Generate the multiple outcomes
exampleDat <- simico_gen_dat(bhFunInv = bhFunInv, obsTimes = 1:3, windowHalf = 0.1, n, p, k, tauSq, gMatCausal, effectSizes)

# Set the initial estimate values
init_beta <-c (rep(c(0, 0, 0, 1, 0), k), 1)

# Run the newton-raphson
nullFit <- simico_fit_null(init_beta = init_beta, epsilon = 10^-5, xDats = exampleDat$fullDat$xDats, lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, k = k, d = d)

# Get the test statistics p-values
out <- simico_out(nullFit = nullFit$beta_fit, xDats = exampleDat$fullDat$xDats, lt_all = exampleDat$leftTimesMat, rt_all = exampleDat$rightTimesMat, Itt = nullFit$jmat, a1 = 1, a2 = 25, G = gMat, k  = k, d = d)

# Print results
# Score statistic
out$multQ

# P-values
out$multP
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
