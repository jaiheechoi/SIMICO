\name{sim_gmat}
\alias{sim_gmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate genetic matrix.
}
\description{
Simulates a n x q genetic matrix.
}
\usage{
sim_gmat(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
% n,q,rho
  \item{n}{Total number of observations.}
  \item{q}{Total number of SNPs.}
  \item{rho}{Common pairwise correlation parameter.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result of a n x q genetic matrix of q SNPs.
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
n = 1500
q = 50
rho = 0.1
gMat <- sim_gmat(n, q, rho)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
