\name{hotperm}
\alias{hotperm}
\alias{print.hotperm}
\alias{summary.hotperm}
\alias{quantile.hotperm}
\alias{print.summary.hotperm}
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
\method{quantile}{hotperm}(x, probs, \dots, lod.thr = NULL)
\method{print}{summary.hotperm}(x, \dots)
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
  \item{x,object}{object of class \code{hotperm} or \code{summary.hotperm}}
  \item{probs}{probability levels for quantiles (\code{1-probs} if all > 0.5);
    default is \code{alpha.levels}}
  \item{lod.thr}{restrict to values above this if not \code{NULL}}
  \item{\dots}{ignored}
}
\author{
  Elias Chaibub Neto and Brian S Yandell
}
\examples{
\dontrun{
set.seed(12345)
hotperm1 <- hotperm(cross = cross1,
                    n.quant = 300,
                    n.perm = 100,
                    lod.thrs = lod.thrs,
                    alpha.levels = alphas,
                    drop.lod = 1.5,
                    verbose = FALSE)
summary(hotperm1)
}
}
\keyword{utilities}
