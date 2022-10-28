\name{ss_fd}
\alias{ss_fd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d/d_{sigma^2}
}
\description{
Calculates the first derivative term with respect to sigma^2.}
\usage{
ss_fd(l, phen, HL_array, HR_array, tpos_all, obs_all, apply_diffs, temp_beta, A_i, no_l_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{Index of first outcome of interest.}
  \item{phen}{List containing all the left and right design matrices.}
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
  \item{tpos_all}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{obs_all}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output is a single value for the first derivative with respect to sigma.}
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
