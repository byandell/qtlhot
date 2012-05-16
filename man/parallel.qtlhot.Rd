\name{parallel.qtlhot}
\alias{parallel.qtlhot}
\alias{qtlhot.phase0}
\alias{big.phase0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Code for parallelizing R/qtlhot.
}
\description{
Code for parallelizing R/qtlhot. See installed parallel directory for
proper use.
}
\usage{
parallel.qtlhot(phase, index = 1, \dots, dirpath = ".")
qtlhot.phase0(dirpath, init.seed = 92387475, len = rep(400,16), n.mar = 185, n.ind = 112,
   n.phe = 100, latent.eff = 0, res.var = 1, lod.thrs, \dots)
big.phase0(dirpath, cross, trait.file, trait.matrix, droptrait.names = NULL,
   keeptrait.names = NULL, lod.thrs, sex = "Sex", trait.index,
   batch.effect = NULL, size.set = 250, offset = 0, subset.sex = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phase}{
%%     ~~Describe \code{phase} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{dirpath}{
%%     ~~Describe \code{dirpath} here~~
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
function (phase, index = 1, ..., dirpath = ".") 
{
    switch(phase, qtlhot.phase1(dirpath, index, ...), qtlhot.phase2(dirpath, 
        index, ...), qtlhot.phase3(dirpath, index, ...), parallel.error(1, 
        phase, index))
    parallel.error(0, phase, index)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
