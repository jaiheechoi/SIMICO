\name{simico_gen_dat}
\alias{simico_gen_dat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
simico_gen_dat()
}
\description{
Generate multiple interval-censored data under proportional hazards model.

}
\usage{
simico_gen_dat(bhFunInv, obsTimes, windowHalf = 0.5, n, k, tauSq, gMatCausal, effectSizes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item {bhFunInv}{the inverse of the baseline hazard function}
  \item{obsTimes}{Vector of the intended observation times}
  \item{windowHalf}{The amount of time before or after the intended obsTimes that a visit might take place}
  \item{n}{number of observations}
  \item{k}{number of outcomes}
  \item{tausq}{variance of subject specific random effect}
  \item{gMatCausal}{matrix of subsetted genetic information for only a select causal SNPs}
  \item{effectSizes}{vector of genetic effects}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  exactTimesMat = exactTimesMat, leftTimesMat = leftTimesMat,
%rightTimesMat = rightTimesMat, obsInd = obsInd, tposInd = tposInd
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
   \item{exactTimesMat}{n x k matrix containing the simulated exact times that the event occurred}
   \item{leftTimesMat}{n x k matrix containing the left time interval that is observed}
   \item{rightTimesMat}{n x k matrix containing the right time interval that is observed}
   \item{obsInd}{n x k matrix containing a indictor for whether that time is right-censored or not}
   \item{tposInd}{n x k matrix containing a indictor for whether that time is left-censored or not}
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
