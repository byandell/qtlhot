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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hot.scan}{
%%     ~~Describe \code{hot.scan} here~~
}
  \item{quant.thr}{
%%     ~~Describe \code{quant.thr} here~~
}
  \item{maps}{
%%     ~~Describe \code{maps} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
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
function (hot.scan, quant.thr = NULL, maps = NULL, main = "") 
{
    for (i in levels(hot.scan$chr)) {
        ylims <- range(c(hot.scan[hot.scan$chr == i, -(1:2)]))
        plot(hot.scan, lodcolumn = 1, chr = i, col = "black", 
            ylab = "hotspot size", ylim = ylims)
        mtext(paste(main, "raw (black) smoothed (red) sliding (blue)"), 
            3, 3)
        plot(hot.scan, lodcolumn = 3, chr = i, col = "red", add = TRUE)
        plot(hot.scan, lodcolumn = 2, chr = i, col = "blue", 
            add = TRUE)
        if (!is.null(maps)) 
            add.rug(i, "", maps, use.cM = TRUE)
        if (!is.null(quant.thr)) 
            abline(h = quant.thr, lwd = 2, lty = 2, col = "red")
    }
    invisible(hot.scan)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
