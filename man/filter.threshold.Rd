\name{filter.threshold}
\alias{filter.threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of threshold results}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filter.threshold(cross, nT, latent.eff, res.var, lod.thrs, drop = 1.5,
Ns, n.perm, alpha.levels, NL.N.thrs, WW.thrs, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
%%     ~~Describe \code{cross} here~~
}
  \item{nT}{
%%     ~~Describe \code{nT} here~~
}
  \item{latent.eff}{
%%     ~~Describe \code{latent.eff} here~~
}
  \item{res.var}{
%%     ~~Describe \code{res.var} here~~
}
  \item{lod.thrs}{
%%     ~~Describe \code{lod.thrs} here~~
}
  \item{drop}{
%%     ~~Describe \code{drop} here~~
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
  \item{NL.N.thrs}{
Results of call to NL.N.permutations.
}
  \item{WW.thrs}{
Results of call to WW.permutations.
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
function (cross, nT, latent.eff, res.var, lod.thrs, drop = 1.5, 
    Ns, n.perm, alpha.levels, verbose = FALSE) 
{
    mycat("scanone", verbose)
    scan <- scanone(cross, pheno.col = c(1:nT), method = "hk")
    mycat("set.to.zero.beyond.drop.int", verbose)
    scan.drop <- set.to.zero.beyond.drop.int(scan, min(lod.thrs), 
        drop)
    mycat("NL.N.permutations", verbose)
    NL.N.thrs <- NL.N.permutations(cross, Ns, n.perm, alpha.levels, 
        lod.thrs, verbose = verbose)
    NL.thrs <- NL.N.thrs[[1]]
    N.thrs <- NL.N.thrs[[2]]
    mycat("WW.permutations", verbose)
    WW.thrs <- WW.permutations(scan.drop, lod.thrs, alpha.levels, 
        n.perm)
    mycat("NL.counts", verbose)
    Nmax <- length(Ns)
    NL <- NL.counts(scan.drop, Nmax, NL.thrs)
    mycat("N.WW.counts", verbose)
    N.WW <- N.WW.counts(scan.drop, lod.thrs, N.thrs, WW.thrs)
    list(NL.thrs = NL.thrs, N.thrs = N.thrs, WW.thrs = WW.thrs, 
        NL = NL, N.counts = N.WW$N, WW.counts = N.WW$WW)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
