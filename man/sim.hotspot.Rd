\name{sim.hotspot}
\alias{sim.hotspot}
\alias{mySimulations}
\alias{sim.null.cross}
\alias{sim.null.pheno.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wrapper routine for simulations.}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.hotspot(nSim, cross, nT, latent.eff, res.var = 1, Ns, n.perm,
  alpha.levels, lod.thrs, drop = 1.5, verbose = FALSE)
mySimulations(\dots)
sim.null.cross(chr.len = rep(400, 16), n.mar = 185, n.ind = 112, type = "bc", n.pheno = 6000, latent.eff = 1.5, res.var = 1, init.seed = 92387475)
sim.null.pheno.data(cross, n.pheno, latent.eff, res.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nSim}{
    Number of simulated sets of phenotypes to create. See details.
}
  \item{cross}{
Object of class \code{cross}. See \code{\link[qtl]{read.cross}}.
}
  \item{nT}{
Number of traits to simulate for cross object.
}
  \item{latent.eff}{
Strength of latent effect, which is included in all traits. See \code{\link{sim.null.cross}}.
}
  \item{res.var}{
Residual variance for traits. Should not affect results.
}
  \item{Ns}{
Vector of sizes of hotspots, typically \code{1:Nmax}, with \code{Nmax}
determined to be large enough to exceed the largest Breitling alpha critical value.
}
  \item{n.perm}{
Number of permutations to perform per realization. Good idea to do 1000,
but this takes time.
}
  \item{alpha.levels}{
Vector of significance levels.
}
  \item{lod.thrs}{
Vector of LOD thresholds, typically single-trait
permutation thresholds for various significance levels.
}
  \item{drop}{
Drop in LOD score examined. LODs below this drop from the maximum for a
chromosome will not be scored.
}
\item{init.seed}{}
\item{len}{}
\item{n.mar}{}
\item{n.ind}{}
\item{n.phe}{}
\item{latent.eff}{}
\item{res.var}{}
\item{lod.thrs}{}
  \item{verbose}{
Verbose output if \code{TRUE}. More detailed output if \code{2}.
}
  \item{\dots}{
Arguments passed directly to \code{sim.hotspot}.
}
}
\details{
Simulate \code{nSim} realizations of cross object with \code{nT} phenotypes with correlation
\code{latent.eff}. All simulations use the same genotypes in the
\code{cross} object.
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
out <- sim.hotspot(100, cross, nT, latent.eff = 1, res.var = 1,
  Ns = 1:700, n.perm = 100,
  alpha.levels = c(1:10)/100, lod.thrs = lod.thrs, drop = 1.5,
  verbose = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
