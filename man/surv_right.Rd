\name{surv_right}
\alias{surv_right}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S_{ik}(R_{ik})
}
\description{
Calculates the survival function of the right time interval for outcome k for quadrature node d.
}
\usage{
surv_right(l, d, temp_beta, phen, r1, k)
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
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output is a vector of the survival values of the left times
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
