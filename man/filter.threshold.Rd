\name{filter.threshold}
\alias{filter.threshold}
\title{Summary of threshold results}
\description{
Summary of threshold results.
}
\usage{
filter.threshold(cross, n.pheno, latent.eff, res.var, lod.thrs, drop.lod = 1.5,
  s.quant, n.perm, alpha.levels, qh.thrs, ww.thrs, verbose = FALSE)
}
\arguments{
  \item{cross}{object of class \code{cross}; see \code{\link[qtl]{read.cross}}}
  \item{n.pheno}{number of traits, or phenotypes}
  \item{latent.eff}{ratio of latent effect SD to residual SD}
  \item{res.var}{residual variance (=SD^2)}
  \item{lod.thrs}{LOD threshold values for range of significance (alpha) levels}
  \item{drop.lod}{LOD drop from max LOD to keep in analysis}
  \item{s.quant}{vector of \code{1:Nmax} with \code{Nmax} the maximum hotspot
  size to be considered}
  \item{n.perm}{number of permutations}
  \item{alpha.levels}{range of significance levels; same length as \code{lod.thrs}}
  \item{qh.thrs}{Results of call to \code{\link{hotperm}}}
  \item{ww.thrs}{Results of call to \code{\link{ww.perm}}}
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
\seealso{\code{\link{hotperm}}, \code{\link{ww.perm}}}
\keyword{utilities}
