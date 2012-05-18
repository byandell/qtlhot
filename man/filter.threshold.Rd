\name{filter.threshold}
\alias{filter.threshold}
\title{Summary of threshold results}
\description{
Summary of threshold results.
}
\usage{
filter.threshold(cross, nT, latent.eff, res.var, lod.thrs, drop = 1.5,
  Ns, n.perm, alpha.levels, NL.N.thrs, WW.thrs, verbose = FALSE)
}
\arguments{
  \item{cross}{object of class \code{cross}; see \code{\link[qtl]{read.cross}}}
  \item{nT}{number of traits}
  \item{latent.eff}{ratio of latent effect SD to residual SD}
  \item{res.var}{residual variance (=SD^2)}
  \item{lod.thrs}{LOD threshold values for range of significance (alpha) levels}
  \item{drop}{LOD drop from max LOD to keep in analysis}
  \item{Ns}{vector of \code{1:Nmax} with \code{Nmax} the maximum hotspot
  size to be considered}
  \item{n.perm}{number of permutations}
  \item{alpha.levels}{range of significance levels; same length as \code{lod.thrs}}
  \item{NL.N.thrs}{Results of call to NL.N.permutations}
  \item{WW.thrs}{Results of call to WW.permutations}
  \item{verbose}{verbose output if \code{TRUE}}
}
\value{
  List with items
  \item{NL.thrs}{}
  \item{N.thrs}{}
  \item{WW.thrs}{}
  \item{NL}{}
  \item{N.counts}{}
  \item{WW.counts}{}
}
\references{
Manichaikul A, Dupuis J, Sen S, Broman KW (2006) Poor performance of
bootstrap confidence intervals for the location of a quantitative trait
locus. Genetics 174: 481-489.
}
\keyword{utilities}
