\name{NL.N.permutations}
\alias{NL.N.permutations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct NL and N permutation tests}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NL.N.permutations(cross, Ns, n.perm, alpha.levels, lod.thrs, drop = 1.5, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
%%     ~~Describe \code{cross} here~~
}
  \item{Ns}{
%%     ~~Describe \code{Ns} here~~
}
  \item{n.perm}{
%%     ~~Describe \code{n.perm} here~~
}
  \item{alpha.levels}{
%%     ~~Describe \code{alpha.levels} here~~
}
  \item{lod.thrs}{
%%     ~~Describe \code{lod.thrs} here~~
}
  \item{drop}{
%%     ~~Describe \code{drop} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (cross, Ns, n.perm, alpha.levels, lod.thrs, drop = 1.5, 
    verbose = FALSE) 
{
    mats <- NL.N.perm(cross, Ns, n.perm, lod.thrs, drop = drop, 
        verbose = verbose)
    NL.N.summary(mats$max.lod.quant, mats$max.N, alpha.levels)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
