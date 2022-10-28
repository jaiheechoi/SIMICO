\name{createInt}
\alias{createInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Time Intervals
}
\description{
Create time intervals around exact observed time.
}
\usage{
createInt(obsTimes, eventTime)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsTimes}{A vector of all the times a subject is observed.}
  \item{eventTime}{The exact event time for the subject.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A 2*1 vector which is the interval of the event time.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Originally from the ICSKAT package.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
obsTimes <- 1:10
eventTime <- 2.5
createInt(obsTimes, eventTime)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
