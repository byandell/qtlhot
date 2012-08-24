\name{hotperm
\alias{hotperm}
\alias{summary.hotperm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Conduct NL and N permutation tests}
\description{
Conduct NL and N permutation tests.
}
\usage{
hotperm(cross, n.quant, n.perm, lod.thrs, alpha.levels, drop.lod = 1.5, verbose = FALSE, init.seed = 0)
\method{print}{hotperm}(x, alpha.levels, \dots)
\method{summary}{hotperm}(x, alpha.levels, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{object of class \code{cross}}
  \item{s.quant}{vector of possible hotspot sizes}
  \item{n.quant}{maximum of \code{s.quant}}
  \item{n.perm}{number of permutations}
  \item{alpha.levels}{vector of significance levels}
  \item{lod.thrs}{vector of LOD thresholds}
  \item{drop.lod}{LOD drop amount for support intervals}
  \item{init.seed}{initial seed for pseudo-random number generation}
  \item{max.lod.quant}{matrix of max LOD (permutation by hotspot size)}
  \item{max.N}{matrix of max count (permutation by LOD threshold)}
  \item{scanmat}{matrix from \code{\link[qtl]{scanone}} without
    \code{chr} and \code{pos} columns}
  \item{x}{object of class \code{ww.perm}}
  \item{chr}{chromosomes identified with rows of \code{scanmat}}
  \item{scan}{object of class \code{\link[qtl]{scanone}}}
  \item{\dots}{ignored}
  \item{verbose}{verbose output if \code{TRUE}}
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
