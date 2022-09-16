\name{fit_null_general}
\alias{fit_null_general}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Null Model
}
\description{
Runs Newton-Raphson to find estimates for beta.
}
\usage{
fit_null_general(init_beta, epsilon, xDats, lt_all, rt_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init_beta}{Starting values of beta for the Newton-Raphson.}
  \item{epsilon}{Stopping criterion for Newton-Raphson.}
  \item{xDats}{List of left and right design matrices and exact times.}
  \item{lt_all}{n x k matrix of left times.}
  \item{rt_all}{n x k matrix of right times}
  \item{k}{Total number of outcomes.}
  \item{d}{Total number of quadrature points.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%temp_beta = temp_beta, iter = iter, diff = diff, jmat = jmat, grad = grad
  \item{temp_beta}{Estimates of the beta coefficients.}
  \item{iter}{Total number of iterations.}
  \item{diff}{Final difference between previous iteration of the beta estimates and the current one.}
  \item{jmat}{I_{\theta\theta} - Information matrix of thetha.}
  \item{grad}{Gradient with respect to theta.}
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
function (init_beta, epsilon, xDats, lt_all, rt_all, k, d)
{
    n = nrow(xDats[[1]]$dmats$right_dmat)
    nocol <- ncol(xDats[[1]]$dmats$right_dmat)
    tpos_all <- matrix(NA, nrow = n, ncol = k)
    obs_all <- matrix(NA, nrow = n, ncol = k)
    for (j in 1:k) {
        tpos_all[, j] <- as.numeric(lt_all[, j] > 0)
        obs_all[, j] <- as.numeric(rt_all[, j] != Inf)
    }
    threedmat <- list()
    for (i in 1:k) {
        threedmat[[i]] <- list(dmats = xDats[[i]]$dmats, lt = lt_all[,
            i], rt = rt_all[, i])
    }
    xAll <- list(xDats = threedmat, ts_all = tpos_all, ob_all = obs_all)
    xDats <- xAll$xDats
    t_all <- xAll$ts_all
    o_all <- xAll$ob_all
    ghDat <- fastGHQuad::gaussHermiteData(d)
    w1 <- ghDat$w
    r1 <- ghDat$x
    iter = 1
    diff = 1
    temp_beta <- matrix(init_beta, nrow = 1)
    while (diff > epsilon & iter < 200) {
        HL_array <- array(c(NA, NA), dim = c(n, d, k))
        HR_array <- array(c(NA, NA), dim = c(n, d, k))
        for (i in 1:k) {
            hl_d <- lapply(1:d, haz_left, l = i, temp_beta = temp_beta,
                phen = xDats[[i]], r1 = r1, k = k)
            hr_d <- lapply(1:d, haz_right, l = i, temp_beta = temp_beta,
                phen = xDats[[i]], r1 = r1, k = k)
            HL_array[, , i] <- simplify2array(hl_d)
            HR_array[, , i] <- simplify2array(hr_d)
        }
        apply_diffs <- array(c(NA, NA), dim = c(n, k, d))
        for (i in 1:k) {
            phen = xDats[[i]]
            apply_diffs[, i, ] <- simplify2array(lapply(1:d,
                surv_diff, l = i, temp_beta = temp_beta, phen = phen,
                r1 = r1, k = k))
        }
        A_i <- get_A(apply_diffs, w1, 100, n)
        A_i[which(A_i == 0)] <- min(A_i[which(A_i != 0)])
        no_l_all <- simplify2array(lapply(1:k, without_one_phen,
            k = k, store = apply_diffs))
        combs <- combn(1:k, 2)
        if (k == 2) {
            no_two_all <- 1
        }
        else {
            no_two_all <- array(data = NA, dim = c(n, d, choose(k,
                2)))
            for (i in 1:choose(k, 2)) {
                no_two_all[, , i] <- without_two_phen(combs[1,
                  i], combs[2, i], k, apply_diffs, n, d)
            }
        }
        grad <- c()
        for (i in 1:k) {
            temp_grad <- fd_term(i, temp_beta, xDats[[i]], d,
                apply_diffs = apply_diffs, A_i = A_i, no_l_all = no_l_all,
                HL_array = HL_array, HR_array = HR_array)
            grad <- c(grad, temp_grad)
        }
        grad_ss <- ss_fd(1, xDats[[1]], HL_array, HR_array, t_all,
            o_all, apply_diffs = apply_diffs, temp_beta = temp_beta,
            A_i = A_i, no_l_all = no_l_all, k = k, d = d)
        grad <- matrix(c(grad, grad_ss), ncol = 1)
        totd <- (nocol * k) + 1
        jmat <- matrix(NA, nrow = totd, ncol = totd)
        for (i in 1:k) {
            d1 <- (nocol * (i - 1)) + 1
            d2 <- nocol * i
            jmat[d1:d2, d1:d2] <- sd_on(i, k, temp_beta, xDats[[i]],
                d, apply_diffs = apply_diffs, A_i = A_i, no_l_all = no_l_all,
                HL_array = HL_array, HR_array = HR_array)
            jmat[d1:d2, totd] <- st_off(i, HL_array, HR_array,
                xAll, apply_diffs = apply_diffs, temp_beta = temp_beta,
                A_i = A_i, no_l_all = no_l_all, no_two_all = no_two_all,
                k = k, d = d)
            jmat[totd, d1:d2] <- t(st_off(i, HL_array, HR_array,
                xAll, apply_diffs = apply_diffs, temp_beta = temp_beta,
                A_i = A_i, no_l_all = no_l_all, no_two_all = no_two_all,
                k = k, d = d))
            idx <- (1:k)[-i]
            for (j in 1:length(idx)) {
                od1 <- (nocol * (idx[j] - 1)) + 1
                od2 <- nocol * idx[j]
                jmat[d1:d2, od1:od2] <- sd_off(i, idx[j], phen_l = xDats[[i]],
                  phen_m = xDats[[idx[j]]], temp_beta, d = d,
                  apply_diffs = apply_diffs, A_i = A_i, HL_array = HL_array,
                  HR_array = HR_array, no_l_all = no_l_all, no_two_all = no_two_all,
                  tpos_all = t_all, obs_all = o_all, k = k)
                jmat[od1:od2, d1:d2] <- sd_off(idx[j], i, phen_l = xDats[[idx[j]]],
                  phen_m = xDats[[i]], temp_beta, d = d, apply_diffs = apply_diffs,
                  A_i = A_i, HL_array = HL_array, HR_array = HR_array,
                  no_l_all = no_l_all, no_two_all = no_two_all,
                  tpos_all = t_all, obs_all = o_all, k = k)
            }
        }
        jmat[totd, totd] <- ss_sd(HL_array, HR_array, xAll, apply_diffs = apply_diffs,
            temp_beta = temp_beta, A_i = A_i, no_l_all = no_l_all,
            no_two_all = no_two_all, k = k, d = d)
        beta_new <- temp_beta - t(grad) \%*\% solve(jmat)
        diff = (-t(grad) \%*\% solve(jmat)) \%*\% t(-t(grad) \%*\%
            solve(jmat))
        iter = iter + 1
        temp_beta <- matrix(beta_new, nrow = 1)
        print(iter)
        print(diff)
    }
    return(list(temp_beta = temp_beta, iter = iter, diff = diff,
        jmat = jmat, grad = grad))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
