\name{pull.highlods}
\alias{pull.highlods}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pull high LOD values with chr and pos.
}
\description{
Pull high LOD values with chr and pos.
}
\usage{
pull.highlods(scans, pheno.col, lod = 4.5, drop.lod = 1.5,
  extend = TRUE, restrict.lod = FALSE)
}
\arguments{
  \item{scans}{object of class \code{\link[qtl]{scanone}}}
  \item{pheno.col}{column numbers for phenotypes (all if missing)}
  \item{lod}{LOD threshold}
  \item{drop.lod}{LOD drop from max to keep for support intervals}
  \item{extend}{extend support interval just past \code{drop.lod};
    matches behavior of \code{\link[qtl]{lodint}} when \code{TRUE}}
  \item{restrict.lod}{restrict to loci above LOD threshold if
    \code{TRUE}, which important when searching for hotspots, but
    probably not when finding peak overlap; matches behavior of
    \code{\link[qtl]{lodint}} when \code{FALSE}}
}
\value{
  Data frame with
  \item{row}{row number in \code{\link[qtl]{scanone}} object}
  \item{phenos}{phenotype column number}
  \item{lod}{LOD score for phenotype at locus indicated by \code{row}} 
}
\author{
Brian S Yandell and Elias Chaibub Neto
}
\examples{
\dontrun{
## Used in Phase2 (big=TRUE) for parallel hotspot runs.
per.scan.hl <- pull.highlods(per.scan, lod = lod.min,
  drop.lod = 1.5, restrict.lod = TRUE)
}
}
\keyword{utilities}
