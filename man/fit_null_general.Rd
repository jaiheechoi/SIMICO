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

TO MAKE
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
