\name{NL.N.permutations}
\alias{NL.N.permutations}
\alias{NL.N.perm}
\alias{NL.N.summary}
\alias{WW.permutations}
\alias{WW.perm}
\alias{WW.summary}
\alias{set.to.zero.beyond.drop.int}
\alias{count.thr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct NL and N permutation tests}
\description{
Conduct NL and N permutation tests.
}
\usage{
NL.N.permutations(cross, Ns, n.perm, alpha.levels, lod.thrs, drop = 1.5, verbose = FALSE)
NL.N.perm(cross, Nmax, n.perm, lod.thrs, drop = 1.5, verbose = FALSE, init.seed = 0)
NL.N.summary(max.lod.quant, max.N, alpha.levels)
WW.permutations(scanmat, lod.thrs, alpha.levels, n.perm, verbose = FALSE)
WW.perm(scanmat, lod.thrs, n.perm, verbose = FALSE)
WW.summary(max.WW, alpha.levels)
set.to.zero.beyond.drop.int(chr, scanmat, thr, drop = 1.5)
count.thr(scan, lod.thrs, droptwo = TRUE)
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
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Elias Chaibub Neto and Brian S Yandell
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
