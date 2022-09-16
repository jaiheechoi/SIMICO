\name{surv_diff}
\alias{surv_diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S_{ik}(L_{ik}) - S_{ik}(R_{ik})
}
\description{
Calculates the difference between the survival functions of the left and right time intervals for outcome k
}
\usage{
surv_diff(l, d, temp_beta, phen, r1, k)
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
The output is a vector of the difference of the survival values of the left times and right times
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
    (surv_left(l, d, temp_beta, phen, r1, k) - surv_right(l,
        d, temp_beta, phen, r1, k))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
