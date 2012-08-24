\name{lod.quantile.perm}
\alias{lod.quantile.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conduct permutations of LOD quantiles.
}
\description{
Conduct permutations of LOD quantiles.
}
\usage{
lod.quantile.permutation(highobj, n.quant, lod.thr, window, chr.pos, n.phe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{highobj}{data frame with \code{row} index and high LOD interval
    retained by \code{\link{pull.highlods}}}
  \item{n.quant}{maximum size of hotspot examined}
  \item{lod.thr}{vector of LOD thresholds}
  \item{window}{size of window (in cM) for running mean}
  \item{chr.pos}{data frame with \code{chr} and \code{pos} including pseudo-markers}
  \item{n.phe}{number of phenotypes}
}
\value{
  List with items
  \item{max.N}{vector of maximum number of traits above each value of \code{lod.thr}}
  \item{max.N.window}{version of \code{max.N} smoothed by \code{window}}
  \item{max.lod.quant}{vector of \code{LOD} quantiles by hotspot size
    from 1 up to \code{min(n.quant, max.N)}}
}
\references{
Chaibub Neto E, Keller MP, Broman AF, Attie AD, Jansen RC, Broman KW,
Yandell BS (2012) Quantile-based permutation thresholds for QTL hotspots.
\emph{Genetics 00}: 000-000.
}
\author{
Brian S Yandell and Elias Chaibub Neto
}
\seealso{\code{\link{hotspot.plot}}, \code{\link{pull.highlods}}
}
\examples{
\dontrun{
lod.sums <- lod.quantile.permutation(scan.hl, Nmax, lod.thr, window,
  chr.pos, n.phe)
}
}
\keyword{utilities}
