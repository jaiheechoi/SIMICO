\name{st_off}
\alias{st_off}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d_{\theta_k}d{\sigma^2}
}
\description{

}
\usage{
st_off(l, HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
  \item{xAll}{List containing the left and right matrices and event times.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{no_two_all}{n x (K - 2) matrix containing the product of apply_diffs across all outcomes K excluding outcomes l and m.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a matrix containing the component of the information matrix of the theta eta parameters for outcome l}
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
function (l, HL_array, HR_array, xAll, apply_diffs, temp_beta,
    A_i, no_l_all, no_two_all, k, d)
{
    ghDat <- fastGHQuad::gaussHermiteData(d)
    w1 <- ghDat$w
    r1 <- ghDat$x
    idx <- (1:k)[-l]
    phen <- xAll$xDats[[l]]
    tpos_all <- xAll$ts_all
    obs_all <- xAll$ob_all
    ldm <- phen$dmats$left_dmat
    rdm <- phen$dmats$right_dmat
    nocol <- ncol(ldm)
    sigmasq <- temp_beta[k * nocol + 1]
    Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)),
        c(1, 3), ifelse(tpos_all == 0, 0, 1), FUN = "*")
    Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)),
        c(1, 3), ifelse(obs_all == 0, 0, 1), FUN = "*")
    Y2[is.na(Y2)] <- 0
    U1 <- sweep((exp(-HL_array) * -HL_array), c(1, 3), ifelse(tpos_all ==
        0, 0, 1), FUN = "*")
    U2 <- sweep((exp(-HR_array) * -HR_array), c(1, 3), ifelse(obs_all ==
        0, 0, 1), FUN = "*")
    U2[is.na(U2)] <- 0
    surv_no_l <- no_l_all
    st_sd <- function(d) {
        if (k == 2) {
            out <- (ldm * (surv_no_l[, d, l] * Y1[, d, l])) -
                (rdm * (surv_no_l[, d, l] * Y2[, d, l])) + ((ldm *
                U1[, d, l]) - (rdm * U2[, d, l])) * (U1 - U2)[,
                d, idx]
        }
        else {
            combs <- combn(1:k, 2)
            with_l_idx <- which(apply(combs, 2, function(x) which(x ==
                l)) == 1 | apply(combs, 2, function(x) which(x ==
                l)) == 2)
            surv_diff_no_two <- no_two_all[, , with_l_idx]
            out <- (ldm * (surv_no_l[, d, l] * Y1[, d, l])) -
                (rdm * (surv_no_l[, d, l] * Y2[, d, l])) + ((ldm *
                U1[, d, l]) - (rdm * U2[, d, l])) * apply((U1 -
                U2)[, d, idx] * surv_diff_no_two[, d, ], 1, sum)
        }
        return(out)
    }
    to_d <- apply(sweep(simplify2array(lapply(1:d, st_sd)), 3,
        w1 * (r1/sqrt(2 * sigmasq)), FUN = "*"), c(1, 2), sum)
    st_t1 <- colSums(sweep(to_d, 1, A_i, FUN = "/"))/sqrt(pi)
    term3_func <- function(d) {
        (ldm * U1[, d, l] * surv_no_l[, d, l]) - (rdm * U2[,
            d, l] * surv_no_l[, d, l])
    }
    t2a <- sweep(apply(sweep(simplify2array(lapply(1:d, term3_func)),
        3, w1, FUN = "*"), c(1, 2), sum), 1, A_i, FUN = "/")/sqrt(pi)
    t2b <- rowSums(sweep(apply((U1 - U2) * surv_no_l/A_i, c(1,
        2), sum), 2, w1 * r1/sqrt(2 * sigmasq), FUN = "*"))/sqrt(pi)
    st_t2 <- colSums(t2a * t2b)
    deriv <- st_t1 - st_t2
    return(deriv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
