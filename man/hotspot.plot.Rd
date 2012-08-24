\name{hotspot.plot}
\alias{qtlhot.scan}
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
qtlhot.scan(cross, scan, max.lod.quant, lod.thrs, probs, level = 0.05)
quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.05)
hotspot.scan(cross, scan.hl, lod.thr, quant.level, window = 5, verbose = FALSE)
hotspot.plot(hot.scan, quant.thr = NULL, maps = NULL, main = "")
sliding.bar.plot(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
}
\arguments{
  \item{max.lod.quant}{matrix of max LOD (permutation by hotspot size)}
  \item{max.N}{matrix of max count (permutation by LOD threshold)}
  \item{max.N.window}{matrix of max count (permutation by LOD threshold)
  for smoothed counts over \code{window}}
  \item{lod.thrs}{vector of LOD thresholds}
  \item{probs}{vector of significance levels for \code{lod.thrs}}
  \item{level}{significance level}
  \item{cross}{object of class \code{cross}}
  \item{scan.hl}{data frame with \code{row} index and high LOD interval
    retained by \code{\link{pull.highlods}}}
  \item{lod.thr}{see \code{lod.thrs}}
  \item{quant.level}{LOD quantiles at significance level?}
  \item{window}{size of window (in cM) for running mean}
  \item{verbose}{verbose output if \code{TRUE}}
  \item{hot.scan}{object returned by \code{hotspot.scan}}
  \item{quant.thr}{quantiles at threshold levels}
  \item{maps}{list of objects of class \code{map}}
  \item{main}{character string as title for plot}
  \item{scan}{object of class \code{scanone}}
  \item{size.thr}{vector of hotspot size thresholds}
  \item{gap}{gap between chromosomes for plot (see \code{\link[qtl]{plot.scanone}})}
  \item{y.axes}{values to plot on Y axis}
}
\value{
The hotspot.scan() result is an object of class \code{scanone}.
}
\seealso{\code{\link{lod.quantile.perm}},
  \code{\link{pull.highlods}}, \code{\link[qtl]{plot.scanone}}} 
\examples{
\dontrun{
## See parallel folder in installed package for details.
## scan.hl is created from multiple permutations.
## max.lod.quant, max.N and max.N.window are elements of lod.sums.
## lod.sums is created by call to lod.quantile.permutation in Phase2 (big=TRUE).
## lod.thrs supplied from earlier calculations; probs correspond to these.
out <- quant.plot(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.05)

hot.scan <- qtlhot.scan(cross, scan, max.lod.quant, lod.thrs, probs, level = 0.05)
plot(hot.scan)
sliding.bar.plot(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
}
}
\keyword{plot}
