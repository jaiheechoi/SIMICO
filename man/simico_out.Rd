\name{simico_out}
\alias{simico_out}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get P-Values
}
\description{
Calculate test statistic and p-values for multiple outcome test and multiple burden test.
}
\usage{
simico_out(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
% nullFit, xDats, lt_all, rt_all, Itt, a1, a2, G, k, d
  \item{nullFit}{Results of the Newton-Raphson: estimates of the beta coefficients.}
  \item{xDats}{List of design matrices.}
  \item{lt_all}{Matrix containing the generated left interval times.}
  \item{rt_all}{Matrix containing the generated right interval times.}
  \item{Itt}{I_{\theta\theta} - Information matrix of thetha.}
  \item{G}{n x q matrix of genetic information.}
  \item{a1}{First shape parameter of beta parameter.}
  \item{a2}{Second shape parameter of beta parameter.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  multQ = gamma_score, multP = pval$Qq, burdQ = burden_score, burdP = p_burden
  \item{multQ}{Score statistic for multiple outcome test.}
  \item{multP}{P-value for multiple outcome test.}
  \item{burdQ}{Score statistic for multiple burden test.}
  \item{burdP}{P-value for multiple burden test.}
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
