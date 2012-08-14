\name{sim.hotspot}
\alias{sim.hotspot}
\alias{mySimulations}
\alias{sim.null.cross}
\alias{sim.null.pheno.data}
\alias{include.hotspots}
\title{Wrapper routine for simulations.}
\description{
Wrapper routine for simulations
}
\usage{
sim.hotspot(nSim, cross, n.phe, latent.eff, res.var = 1, Ns, n.perm,
  alpha.levels, lod.thrs, drop = 1.5, verbose = FALSE)
mySimulations(\dots)
sim.null.cross(chr.len = rep(400, 16), n.mar = 185, n.ind = 112,
  type = "bc", n.phe = 6000, latent.eff = 1.5, res.var = 1,
  init.seed = 92387475)
sim.null.pheno.data(cross, n.phe, latent.eff, res.var)
include.hotspots(cross, hchr, hpos, hsize, Q.eff, latent.eff,
  lod.range.1, lod.range.2, lod.range.3, res.var=1, n.phe, init.seed)
}
\arguments{
  \item{nSim}{Number of simulated sets of phenotypes to create. See details.}
  \item{cross}{Object of class \code{cross}. See \code{\link[qtl]{read.cross}}.}
  \item{n.phe}{Number of traits, or phenotypes, to simulate for cross object.}
  \item{latent.eff}{Strength of latent effect, which is included in all traits. See \code{\link{sim.null.cross}}.}
  \item{res.var}{Residual variance for traits. Should not affect results.}
  \item{Ns}{Vector of sizes of hotspots, typically \code{1:Nmax}, with \code{Nmax}
    determined to be large enough to exceed the largest Breitling alpha critical value.}
  \item{n.perm}{Number of permutations to perform per realization. Good idea to do 1000,
    but this takes time.}
  \item{alpha.levels}{Vector of significance levels.}
  \item{lod.thrs}{Vector of LOD thresholds, typically single-trait
    permutation thresholds for various significance levels.}
  \item{drop}{Drop in LOD score examined. LODs below this drop from the maximum for a
    chromosome will not be scored.}
  \item{init.seed}{initial seed for pseudo-random number generation}
  \item{chr.len}{vector of chromosome lengths}
  \item{n.mar}{number of markers}
  \item{n.ind}{number of individuals}
  \item{type}{type of cross}
  \item{hchr,hpos,hsize}{vectors for hotspot chromosomes, positions, and
  sizes}
  \item{Q.eff}{QTL effect}
  \item{lod.range.1,lod.range.2,lod.range.3}{2-vectors of LOD ranges for
  multiple purposes}
  \item{verbose}{Verbose output if \code{TRUE}. More detailed output if \code{2}.}
  \item{\dots}{Arguments passed directly to \code{sim.hotspot}.}
}
\details{
Simulate \code{nSim} realizations of cross object with \code{n.phe} phenotypes with correlation
\code{latent.eff}. All simulations use the same genotypes in the
\code{cross} object.
}
\value{
  \code{sim.null.cross} simulates an object of class \code{cross}.
  \code{sim.null.pheno.data} simulates a data frame of phenotypes.
  \code{sim.hotspot} uses these other routines to simulate a hotspot,
  returning an list object.
}
\author{
Elias Chaibub Neto and Brian S. Yandell
}
\seealso{
\code{\link{sim.null.cross}}, \code{\link[qtl]{read.cross}}.
}
\examples{
\dontrun{
## This will take DAYS to run.
cross <- sim.null.cross()
lod.thrs <- c(4.63,4.17,3.93,3.76,3.65,3.56,3.47,3.39,3.34,3.30)
out <- sim.hotspot(100, cross, n.phe, latent.eff = 1, res.var = 1,
  Ns = 1:700, n.perm = 100,
  alpha.levels = c(1:10)/100, lod.thrs = lod.thrs, drop = 1.5,
  verbose = TRUE)
}
}
\keyword{utilities}
