\name{get_A}
\alias{get_A}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get A vector
}
\description{
Product of difference of survival terms of the left and right interval times, across all outcomes K, summed over all quadrature nodes D
}
\usage{
get_A(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
% store, weights, d, n
  \item{store}{Matrix of difference of survival values of the left and right time intervals.}
  \item{weights}{Gaussian quadrature weights.}
  \item{d}{Total number of Gaussian quadrature nodes.}
  \item{n}{Total number of observations.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a vector used to compute the derivative terms
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
\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
