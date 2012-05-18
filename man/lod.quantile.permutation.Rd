\name{lod.quantile.permutation}
\alias{lod.quantile.permutation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conduct permutations of LOD quantiles.
}
\description{
Conduct permutations of LOD quantiles.
}
\usage{
lod.quantile.permutation(cat.scan.hl, N, lod.thr, window, chr.pos, n.phe)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cat.scan.hl}{}
  \item{N}{}
  \item{lod.thr}{}
  \item{window}{}
  \item{chr.pos}{}
  \item{n.phe}{}
}
\value{
  List with items
  \item{max.lod.quant}{}
  \item{max.N}{list containing \code{max.N} and \code{max.N.window}}
}
\references{
Chaibub Neto E, Keller MP, Broman AF, Attie AD, Jansen RC, Broman KW,
Yandell BS (2012) Quantile-based permutation thresholds for QTL hotspots.
\emph{Genetics 00}: 000-000.
}
\author{
Brian S Yandell and Elias Chaibub Neto
}
\seealso{\code{\link{hotspot.plot}}
}
\examples{
\dontrun{
lod.sums <- lod.quantile.permutation(scan.hl, Nmax, lod.thr, window,
  chr.pos, n.phe)
}
}
\keyword{utilities}
