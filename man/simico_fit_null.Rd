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
  \item{init_beta}{starting values for NR}
  \item{epsilon}{stopping criterion for NR}
  \item{xDats}{list of design matrices}
  \item{lt_all}{n x k matrix of left times}
  \item{rt_all}{n x k matrix of right times}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  temp_beta = temp_beta, iter = iter, diff = diff, jmat = jmat, grad = grad
 \item{temp_beta}{vector of fitted coefficients}
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
ADD LATER
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
