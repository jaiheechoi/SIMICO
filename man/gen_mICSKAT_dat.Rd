\name{gen_mICSKAT_dat}
\alias{gen_mICSKAT_dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
gen_mICSKAT_dat
}
\description{
Generate multiple interval-censored data under proportional hazards model.
}
\usage{
gen_mICSKAT_dat(bhFunInv, obsTimes, windowHalf = 0.5, n, k, tauSq, gMatCausal, effectSizes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item {bhFunInv}{The inverse of the baseline hazard function.}
  \item{obsTimes}{Vector of the intended observation times.}
  \item{windowHalf}{The amount of time before or after the intended obsTimes that a visit might take place.}
  \item{n}{Number of observations.}
  \item{k}{Number of outcomes.}
  \item{tausq}{Variance of subject specific random effect.}
  \item{gMatCausal}{Matrix of subsetted genetic information for only a select causal SNPs.}
  \item{effectSizes}{Vector of genetic effects.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
% exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
%              rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd, xAll = xAll
  \item{exactTimesMat}{Matrix containing the generated exact times.}
  \item{leftTimesMat}{Matrix containing the generated left interval times.}
  \item{rightTimesMat}{Matrix containing the generated right interval times.}
  \item{obsInd}{n x k matrix containing a indictor for whether that time is right-censored or not.}
  \item{tposInd}{n x k matrix containing a indictor for whether that time is left-censored or not.}
  \item{xAll}{List containing all the design matrices, left and right times, and the censoring terms.}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
