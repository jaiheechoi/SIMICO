\name{sd_on}
\alias{sd_on}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d{\sigma^2}^2
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sd_on(l, k, temp_beta, phen, d, apply_diffs, A_i, no_l_all, HL_array, HR_array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{phen}{List containing the left and right design matrices and interval times for outcome l.}
  \item{d}{Total number of quadrature nodes.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
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
function (l, k, temp_beta, phen, d, apply_diffs, A_i, no_l_all,
    HL_array, HR_array)
{
    ghDat <- fastGHQuad::gaussHermiteData(d)
    w1 <- ghDat$w
    r1 <- ghDat$x
    no_l <- no_l_all[, , l]
    left_dmat <- phen$dmats$left_dmat
    right_dmat <- phen$dmats$right_dmat
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)
    hl_d <- HL_array[, , l]
    hr_d <- HR_array[, , l]
    get_sd <- function(hl_d, hr_d, phen, d, no_l) {
        left_dmat <- phen$dmats$left_dmat
        right_dmat <- phen$dmats$right_dmat
        lt <- phen$lt
        rt <- phen$rt
        tpos_ind <- as.numeric(lt > 0)
        obs_ind <- as.numeric(rt != Inf)
        ul_1 <- ifelse(tpos_ind == 0, 0, -hl_d[, d] * exp(-hl_d[,
            d]) + (hl_d[, d]^2 * exp(-hl_d[, d])))
        ur_1 <- ifelse(obs_ind == 0, 0, -hr_d[, d] * exp(-hr_d[,
            d]) + (hr_d[, d]^2 * exp(-hr_d[, d])))
        ur_1[which(is.na(ur_1))] <- 0
        sd_term1 <- t(left_dmat) \%*\% ((no_l[, d] * as.numeric(ul_1/A_i)) *
            left_dmat)
        sd_term2 <- t(right_dmat) \%*\% ((no_l[, d] * as.numeric(ur_1)/A_i) *
            right_dmat)
        sd_5x5 <- (sd_term1 - sd_term2) * w1[d]
        return(sd_5x5)
    }
    derivs <- lapply(1:d, get_sd, hl_d = hl_d, hr_d = hr_d, phen = phen,
        no_l = no_l)
    derivs_array <- simplify2array(derivs)
    term1 <- apply(derivs_array, c(1, 2), sum)/sqrt(pi)
    first_deriv <- function(l, d, hl_d, hr_d, phen) {
        left_dmat <- phen$dmats$left_dmat
        right_dmat <- phen$dmats$right_dmat
        lt <- phen$lt
        rt <- phen$rt
        tpos_ind <- as.numeric(lt > 0)
        obs_ind <- as.numeric(rt != Inf)
        U1 <- left_dmat * ifelse(tpos_ind == 0, 0, (exp(-hl_d[,
            d]) * -hl_d[, d]/A_i))
        U2 <- right_dmat * ifelse(obs_ind == 0, 0, (exp(-hr_d[,
            d]) * -hr_d[, d]/A_i))
        U2[is.na(U2)] <- 0
        inside <- U1 - U2
        return(inside)
    }
    insides <- lapply(1:d, first_deriv, l = l, hl_d = hl_d, hr_d = hr_d,
        phen = phen)
    mult_together <- function(d, arrayA, listB, weights) {
        arrayA[, d] * listB[[d]] * weights[d]
    }
    deriv_prod <- lapply(1:d, mult_together, arrayA = no_l, listB = insides,
        weights = w1)
    dp_array <- simplify2array(deriv_prod)
    term2_a <- apply(dp_array, c(1, 2), sum)/sqrt(pi)
    term2 <- (t(term2_a) \%*\% term2_a)
    sd_on <- term1 - term2
    return(sd_on)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
