\name{without_one_phen}
\alias{without_one_phen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival Difference Product without Outcome l
}
\description{
Calculate the product of the difference between survival terms excluding that of the outcome of interest.
}
\usage{
without_one_phen(l, k, store)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{store}{Array of difference between left and right survival values.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A n x (k-1) matrix where each column is the product of all the differences of left and right survival values across all outcomes excluding the column index outcome.
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
