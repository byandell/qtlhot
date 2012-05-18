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
  \item{cross}{}
  \item{Ns}{}
  \item{Nmax}{}
  \item{n.perm}{}
  \item{alpha.levels}{}
  \item{lod.thrs}{}
  \item{drop}{}
  \item{init.seed}{}
  \item{max.lod.quant}{}
  \item{max.N}{}
  \item{scanmat}{}
  \item{max.WW}{}
  \item{chr}{}
  \item{thr}{}
  \item{droptwo}{}
  \item{scan}{}
  \item{verbose}{verbose output if \code{TRUE}}
}
\author{
Elias Chaibub Neto and Brian S Yandell
}
\examples{
}
\keyword{utilities}
