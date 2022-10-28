\name{gamma_off}
\alias{gamma_off}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d_{gamma_k}d{gamma_j}
}
\description{
Calculates the [off-diagonal] Information matrix term for I_{gamma gamma} with respect to outcome k and outcome j
}
\usage{
gamma_off(l, m, HL_array, HR_array, tpos_all, obs_all, temp_beta, A_i, no_l_all, no_two_all, gMat, a1, a2, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest}
  \item{m}{Index of second outcome of interest}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not}
  \item{temp_beta}{vector of fitted coefficients}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l}
  \item{no_two_all}{n x (K - 2) matrix containing the product of apply_diffs across all outcomes K excluding the outcomes l and m}
  \item{gMat}{n x q matrix of genetic information}
  \item{a1}{First shape parameter of beta parameter}
  \item{a2}{Second shape parameter of beta parameter}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a matrix containing the component of the information matrix of the gamma parameter for outcomes l and m
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
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
