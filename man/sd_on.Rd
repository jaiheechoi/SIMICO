\name{sd_on}
\alias{sd_on}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d{sigma^2}^2
}
\description{
Calculates the Information matrix term of I_{sigma^2 sigma^2} for outcome l.
}
\usage{
sd_on(l, k, temp_beta, phen, d, apply_diffs, A_i,
   no_l_all, HL_array, HR_array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{k}{Total number of outcomes.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{phen}{List containing the left and right design matrices and interval times for outcome l.}
  \item{d}{Total number of quadrature nodes.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output is a single value for the second derivative with respect to sigma.}
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
