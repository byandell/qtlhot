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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scans}{
%%     ~~Describe \code{scans} here~~
}
  \item{pheno.col}{
%%     ~~Describe \code{pheno.col} here~~
}
  \item{lod}{
%%     ~~Describe \code{lod} here~~
}
  \item{drop.lod}{
%%     ~~Describe \code{drop.lod} here~~
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
function (scans, pheno.col, lod = 4.5, drop.lod = 1.5) 
{
    if (missing(pheno.col)) 
        pheno.col <- names(scans)[-(1:2)]
    x <- as.matrix(scans[, -(1:2), drop = FALSE])
    wh <- which(x >= lod)
    rr <- row(x)[wh]
    cc <- col(x)[wh]
    lod <- x[wh]
    chr <- scans$chr[rr]
    tmp <- interaction(cc, chr, drop = TRUE)
    maxlod <- tapply(lod, tmp, max)
    wh <- which(maxlod[tmp] <= lod + drop.lod)
    rr <- rr[wh]
    cc <- cc[wh]
    lod <- lod[wh]
    cbind.data.frame(row = rr, phenos = pheno.col[cc], lod = lod)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
