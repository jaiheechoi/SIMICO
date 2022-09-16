\name{gamma_on}
\alias{gamma_on}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d_{gamma_k}d{gamma_k}
}
\description{
Calculates the [on-diagonal] Information matrix term for I_{\gamma\gamma} with respect to outcome k
}
\usage{
gamma_on(l, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, gMat, a1, a2, d)
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
The output is a matrix containing the component of the information matrix of the gamma parameter for outcome l
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
    Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)),
        c(1, 3), ifelse(tpos_all == 0, 0, 1), FUN = "*")
    Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)),
        c(1, 3), ifelse(obs_all == 0, 0, 1), FUN = "*")
    Y2[is.na(Y2)] <- 0
    term1 <- function(d, l) {
        (t(Z_w * (Y1 - Y2)[, d, l] * no_l_all[, d, l]/A_i)) \%*\%
            Z_w
    }
    gt_t1 <- apply(sweep(simplify2array(lapply(1:d, term1, l = l)),
        3, w1, FUN = "*"), c(1, 2), sum)/sqrt(pi)
    U1 <- sweep((exp(-HL_array) * -HL_array), c(1, 3), ifelse(tpos_all ==
        0, 0, 1), FUN = "*")
    U2 <- sweep((exp(-HR_array) * -HR_array), c(1, 3), ifelse(obs_all ==
        0, 0, 1), FUN = "*")
    U2[is.na(U2)] <- 0
    no_l <- no_l_all[, , l]
    fd_t <- function(l, d) {
        out <- sweep(Z_w, 1, (U1 - U2)[, d, l] * no_l[, d]/A_i,
            FUN = "*")
        return(out)
    }
    to_d <- apply(sweep(simplify2array(lapply(1:d, fd_t, l = l)),
        3, w1, FUN = "*"), c(1, 2), sum)/sqrt(pi)
    gt_t2 <- t(t(to_d) \%*\% to_d)
    return(gt_t1 - gt_t2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
