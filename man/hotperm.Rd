\name{hotperm}
\alias{hotperm}
\alias{print.hotperm}
\alias{summary.hotperm}
\alias{quantile.hotperm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct NL and N permutation tests}
\description{
Conduct NL and N permutation tests.
}
\usage{
hotperm(cross, n.quant, n.perm, lod.thrs, alpha.levels, drop.lod = 1.5,
  window = NULL, verbose = FALSE, init.seed = 0)
\method{print}{hotperm}(x, \dots)
\method{summary}{hotperm}(object, \dots)
\method{quantile}{hotperm}(x, probs = 0.95, \dots, lod.thr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{object of class \code{cross}}
  \item{n.quant}{maximum of \code{s.quant}}
  \item{n.perm}{number of permutations}
  \item{lod.thrs}{vector of LOD thresholds}
  \item{alpha.levels}{vector of significance levels}
  \item{drop.lod}{LOD drop amount for support intervals}
  \item{window}{window size for smoothed hotspot size}
  \item{verbose}{verbose output if \code{TRUE}}
  \item{init.seed}{initial seed for pseudo-random number generation}
  \item{x,object}{object of class \code{hotperm}}
  \item{probs}{probability levels for quantiles (should be > 0.5)}
  \item{lod.thr}{restrict to values above this if not \code{NULL}}
  \item{\dots}{ignored}
}
\author{
Elias Chaibub Neto and Brian S Yandell
}
\examples{
\dontrun{
qh <- hotperm(cross, n.quant, n.perm, lod.thrs, drop.lod = 1.5,
  verbose = FALSE, init.seed = 0)
summary(qh)
}
}
\keyword{utilities}
