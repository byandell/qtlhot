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
  \item{cat.scan.hl}{
%%     ~~Describe \code{cat.scan.hl} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{lod.thr}{
%%     ~~Describe \code{lod.thr} here~~
}
  \item{window}{
%%     ~~Describe \code{window} here~~
}
  \item{chr.pos}{
%%     ~~Describe \code{chr.pos} here~~
}
  \item{n.phe}{
%%     ~~Describe \code{n.phe} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cat.scan.hl, N, lod.thr, window, chr.pos, n.phe) 
{
    n.chr <- levels(chr.pos$chr)
    quant <- get.tails(cat.scan.hl, n.quant = N)
    N <- ncol(quant)
    max.lod.quant <- apply(quant, 2, max, na.rm = T)
    names(max.lod.quant) <- paste(paste(round(1 - as.numeric(dimnames(quant)[[2]])/n.phe, 
        4) * 100, "\%", sep = ""), 1:N, sep = "_")
    max.N <- make.max.N(cat.scan.hl, lod.thr, chr.pos, window)
    list(max.lod.quant = max.lod.quant, max.N = max.N)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
