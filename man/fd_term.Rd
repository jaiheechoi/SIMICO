\name{fd_term}
\alias{fd_term}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d/d\theta_l
}
\description{
Calculate the first derivative of the theta terms for outcome l
}
\usage{
fd_term(l, temp_beta, phen,d, apply_diffs, A_i, no_l_all,HL_array, HR_array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%l, temp_beta, phen,d, apply_diffs, A_i, no_l_all,HL_array, HR_array
  \item{l}{Outcome of interest.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{phen}{list containing the covariate design matrices.}
  \item{d}{Number of quadrature nodes.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  fd
The output is a 1 x (p + 2) vector of the first derivative terms for outcome l.
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
