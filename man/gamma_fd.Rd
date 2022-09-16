\name{gamma_fd}
\alias{gamma_fd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d/d_{gamma_k}
\description{
Calculates the gradient term for U_g for the score statistic.
}
\usage{
gamma_fd(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{gMat}{n x q matrix of genetic information}
  \item{a1}{First shape parameter of beta parameter}
  \item{a2}{Second shape parameter of beta parameter}
  \item{d}{Number of quadrature nodes}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a vector containing the first derivative with respect to gamma.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (l, HL_array, HR_array, tpos_all, obs_all, temp_beta,
    A_i, no_l_all, gMat, a1, a2, d)
{
    ghDat <- fastGHQuad::gaussHermiteData(d)
    w1 <- ghDat$w
    r1 <- ghDat$x
    MAF <- apply(gMat, 2, function(x) mean(x, na.rm = T)/2)
    beta <- dbeta(MAF, a1, a2, ncp = 0, log = FALSE)
    Z_w = t(t(gMat) * (beta))
    U1 <- sweep((exp(-HL_array) * -HL_array), c(1, 3), ifelse(tpos_all ==
        0, 0, 1), FUN = "*")
    U2 <- sweep((exp(-HR_array) * -HR_array), c(1, 3), ifelse(obs_all ==
        0, 0, 1), FUN = "*")
    U2[is.na(U2)] <- 0
    no_l <- no_l_all[, , l]
    t1 <- t(Z_w) \%*\% ((U1 - U2)[, , l] * no_l/A_i)
    deriv <- apply(sweep(t1, 2, w1, FUN = "*"), 1, sum)/sqrt(pi)
    return(deriv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
