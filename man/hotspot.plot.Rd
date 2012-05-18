\name{hotspot.plot}
\alias{hotspot.scan}
\alias{hotspot.plot}
\alias{quant.plot}
\alias{sliding.bar.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot routines for hotspots.
}
\description{
Plot hotspots across genome.
}
\usage{
quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)
hotspot.scan(cross, scan.hl, lod.thr, quant.level, window = 5, verbose = FALSE)
hotspot.plot(hot.scan, quant.thr = NULL, maps = NULL, main = "")
sliding.bar.plot(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
}
\arguments{
  \item{max.lod.quant}{}
  \item{max.N}{}
  \item{max.N.window}{}
  \item{lod.thrs}{}
  \item{probs}{}
  \item{level}{}
  \item{cross}{}
  \item{scan.hl}{}
  \item{lod.thr}{}
  \item{quant.level}{}
  \item{window}{}
  \item{verbose}{verbose output if \code{TRUE}}
  \item{hot.scan}{}
  \item{quant.thr}{}
  \item{maps}{}
  \item{main}{}
  \item{scan}{}
  \item{size.thr}{}
  \item{gap}{}
  \item{y.axes}{}
}
\value{
The hotspot.scan() result is an object of class \code{scanone}.
}
\seealso{\code{\link{lod.quantile.permutation}}, \code{\link[qtl]{plot.scanone}}}
\examples{
\dontrun{
## See parallel folder in installed package for details.
## scan.hl is created from multiple permutations.
## max.lod.quant, max.N and max.N.window are elements of lod.sums.
## lod.sums is created by call to lod.quantile.permutation in Phase2 (big=TRUE).
## lod.thrs supplied from earlier calculations; probs correspond to these.
load("Phase3.RData")
out <- quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)

hot.scan <- hotspot.scan(cross, scan.hl, out$lod.thr, out$quant.level)
hotspot.plot(hot.scan)
sliding.bar.plot(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
}
}
\keyword{plot}
