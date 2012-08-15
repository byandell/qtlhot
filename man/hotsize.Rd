\name{hotsize}
\alias{hotsize}
\alias{hotsize.scanone}
\alias{hotsize.highlod}
\alias{print.hotsize}
\alias{summary.hotsize}
\alias{plot.hotsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hotspot size routines.
}
\description{
Determine hotspot sizes and display. Use individual threshold and
quantile thresholds as provided.
}
\usage{
hotsize(hotobject, \dots)
\method{hotsize}{scanone}(hotobject, lod.thr = NULL, drop.lod = 1.5, \dots)
\method{hotsize}{highlod}(hotobject, lod.thr = NULL, window = NULL,
  quant.level = NULL, \dots)
\method{print}{hotsize}(object, \dots)
\method{summary}{hotsize}(object, \dots)
\method{plot}{hotsize}(x, \dots)
}
\arguments{
  \item{hotobject}{object of class \code{\link[qtl]{scanone}} or \code{\link{highlod}}}
  \item{lod.thr}{LOD threshold}
  \item{drop.lod}{LOD drop from max to keep for support intervals}
  \item{window}{window width in cM for smoothing hotspot size; not used
    if \code{0} or \code{NULL}}
  \item{quant.level}{vector of LOD levels for 1 up to
    \code{length(quant.level)} size hotspots}
  \item{x,object}{object of class \code{hotsize}}
  \item{\dots}{arguments passed along}
}
\value{
  \code{hotsize} methods return an object of class \code{hotsize}, which
  is essentially an object of class \code{\link[qtl]{summary.scanone}}
  with additional attributes for \code{lod.thr}, \code{window}, and
  \code{quant.level}. 
}
\author{
Brian S Yandell and Elias Chaibub Neto
}
\examples{
\dontrun{
scan1 <- scanone(cross1, pheno.col = 1:1000, method = "hk")
hots1 <- hotsize(scan1, lod.thr = 3)
summary(hots1)
}
}
\keyword{utilities}
