\name{parallel.qtlhot}
\alias{parallel.qtlhot}
\alias{qtlhot.phase0}
\alias{big.phase0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Code for parallelizing R/qtlhot.
}
\description{
Code for parallelizing R/qtlhot. See installed parallel directory for
proper use.
}
\usage{
parallel.qtlhot(phase, index = 1, \dots, dirpath = ".")
qtlhot.phase0(dirpath, init.seed = 92387475, len = rep(400,16), n.mar = 185, n.ind = 112,
   n.phe = 100, latent.eff = 0, res.var = 1, lod.thrs, \dots)
big.phase0(dirpath, cross, trait.file, trait.matrix, droptrait.names = NULL,
   keeptrait.names = NULL, lod.thrs, sex = "Sex", trait.index,
   batch.effect = NULL, size.set = 250, offset = 0, subset.sex = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phase}{}
  \item{index}{}
  \item{\dots}{}
  \item{dirpath}{}
  \item{init.seed}{}
  \item{len}{}
  \item{n.mar}{}
  \item{n.ind}{}
  \item{n.phe}{}
  \item{latent.eff}{}
  \item{res.var}{}
  \item{lod.thrs}{}
  \item{cross}{}
  \item{trait.file}{}
  \item{trait.matrix}{}
  \item{droptrait.names}{}
  \item{keeptrait.names}{}
  \item{sex}{}
  \item{trait.index}{}
  \item{batch.effect}{}
  \item{size.set}{}
  \item{offset}{}
  \item{subset.sex}{}
  \item{verbose}{verbose output if \code{TRUE}}
}
\author{
Brian S Yandell and Elias Chaibub Neto
}
\keyword{utilities}
