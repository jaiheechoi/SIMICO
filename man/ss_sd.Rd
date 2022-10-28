\name{ss_sd}
\alias{ss_sd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
d^2/d_{sigma^2}^2
}
\description{
Calculates the second derivative term with respect to sigma^2.}
\usage{
ss_sd(HL_array, HR_array, xAll, apply_diffs, temp_beta, A_i, no_l_all, no_two_all, k, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HL_array}{n x K matrix containing all the hazard values for the left times.}
  \item{HR_array}{n x K matrix containing all the hazard values for the right times.}
  \item{xAll}{List containing the left and right matrices and event times.}
  \item{apply_diffs}{Matrix containing the differences between survival functions of the left and right time intervals.}
  \item{temp_beta}{vector of fitted coefficients.}
  \item{A_i}{Product of apply_diffs across all outcomes K summed over all quadrature nodes D.}
  \item{no_l_all}{n x (K - 1) matrix containing the product of apply_diffs across all outcomes K excluding the current outcome l.}
  \item{no_two_all}{n x (K - 2) matrix containing the product of apply_diffs across all outcomes K excluding outcomes l and m.}
  \item{k}{Total number of outcomes.}
  \item{d}{Number of quadrature nodes.}

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
