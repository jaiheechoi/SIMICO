\name{without_two_phen}
\alias{without_two_phen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival Difference Product without Outcomes l and m
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
without_two_phen(l, m, k, store, n, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{The first outcome of interest.}
  \item{m}{The second outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{store}{Array of difference between left and right survival values.}
  \item{n}{Total number of observation.}
  \item{d}{Total number of quadrature nodes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A n x (k-2) matrix containing the product of all the differences of left and right survival values across all outcomes excluding outcomes l and m.
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
