\name{ss_sd}
\alias{ss_sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d_{\sigma^2}^2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ss_sd(HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
The output is a single value for the second derivative with respect to sigma.}
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
function (HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i,
    no_l_all, no_two_all, k, d)
{
    ghDat <- fastGHQuad::gaussHermiteData(d)
    w1 <- ghDat$w
    r1 <- ghDat$x
    nocol <- ncol(xAll$xDats[[1]]$dmats$right_dmat)
    sigmasq <- temp_beta[k * nocol + 1]
    tpos_all <- xAll$ts_all
    obs_all <- xAll$ob_all
    U1 <- sweep((exp(-HL_array) * -HL_array), c(1, 3), ifelse(tpos_all ==
        0, 0, 1), FUN = "*")
    U2 <- sweep((exp(-HR_array) * -HR_array), c(1, 3), ifelse(obs_all ==
        0, 0, 1), FUN = "*")
    U2[is.na(U2)] <- 0
    term_3_rep <- function(l) {
        idx <- (1:k)[-l]
        combs <- combn(1:k, 2)
        to_idx <- function(m) {
            min_k <- min(l, m)
            max_k <- max(l, m)
            if (k == 2) {
                out <- (U1 - U2)[, , l] * (U1 - U2)[, , m]
            }
            else {
                out <- (U1 - U2)[, , l] * (U1 - U2)[, , m] *
                  no_two_all[, , which(combs[1, ] == min_k &
                    combs[2, ] == max_k)]
            }
            return(out)
        }
        apply_idx <- simplify2array(lapply(idx, to_idx))
        return(apply(apply_idx, c(1, 2), sum))
    }
    term_three <- simplify2array(lapply(1:k, term_3_rep))
    surv_no_l <- no_l_all
    Y1 <- sweep(-HL_array * exp(-HL_array) + (HL_array^2 * exp(-HL_array)),
        c(1, 3), ifelse(tpos_all == 0, 0, 1), FUN = "*")
    Y2 <- sweep(-HR_array * exp(-HR_array) + (HR_array^2 * exp(-HR_array)),
        c(1, 3), ifelse(obs_all == 0, 0, 1), FUN = "*")
    Y2[is.na(Y2)] <- 0
    S1 <- sweep((exp(-HL_array) * -HL_array), c(1, 3), ifelse(tpos_all ==
        0, 0, 1), FUN = "*")
    S2 <- sweep((exp(-HR_array) * -HR_array), c(1, 3), ifelse(obs_all ==
        0, 0, 1), FUN = "*")
    S2[is.na(S2)] <- 0
    score <- sweep((U1 - U2) * surv_no_l, 2, -r1/(2 * sqrt(2 *
        sigmasq)), FUN = "*") + sweep((Y1 - Y2) * surv_no_l,
        2, (r1/sqrt(2 * sigmasq))^2, FUN = "*") + sweep(term_three,
        2, (r1/sqrt(2 * sigmasq))^2, FUN = "*")
    out_sumk <- apply(score, c(1, 2), sum)
    sum_d <- apply((t(out_sumk) * w1), 2, sum)
    ss_sd_t1 <- sum(sum_d/A_i)/sqrt(pi)
    ss_sd_t2 <- sum((apply(t(apply((S1 - S2) * surv_no_l, c(1,
        2), sum)) * w1 * r1/sqrt(2 * sigmasq), 2, sum)/A_i/sqrt(pi))^2)
    return(ss_sd_t1 - ss_sd_t2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
