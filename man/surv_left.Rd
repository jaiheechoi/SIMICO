\name{surv_left}
\alias{surv_left}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
surv_left(l, d, temp_beta, phen, r1, k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Outcome of interest.}
  \item{d}{Total number of Gaussian quadrature nodes.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{phen}{list of data matrices containing both left and right information.}
  \item{r1}{Gaussian quadrature nodes.}
  \item{k}{Total number of outcomes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a vector of the survival values of the left times
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
function (l, d, temp_beta, phen, r1, k)
{
    covcol <- ncol(phen$dmats$right_dmat)
    sub_beta <- temp_beta[((l - 1) * covcol + 1):(l * covcol)]
    sigmasq <- temp_beta[k * covcol + 1]
    lt <- phen$lt
    rt <- phen$rt
    tpos_ind <- as.numeric(lt > 0)
    obs_ind <- as.numeric(rt != Inf)
    left_dmat <- phen$dmats$left_dmat
    hl1 <- as.numeric(exp(left_dmat \%*\% t(matrix(sub_beta, nrow = 1)) +
        sqrt(2 * sigmasq) * r1[d]))
    sl1 <- ifelse(tpos_ind == 0, 1, exp(-hl1))
    return(sl1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
