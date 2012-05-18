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
pull.highlods(scans, pheno.col, lod = 4.5, drop.lod = 1.5)
}
\arguments{
  \item{scans}{}
  \item{pheno.col}{}
  \item{lod}{}
  \item{drop.lod}{}
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
## Used in Phase2 (big=TRUE) for parallel runs.
per.scan.hl <- pull.highlods(per.scan, lod = lod.min, drop.lod = drop.lod)
}
}
\keyword{utilities}
