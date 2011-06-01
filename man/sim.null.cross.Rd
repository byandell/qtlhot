\name{sim.null.cross}
\alias{sim.null.cross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate cross with no correlation}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.null.cross(chr.len = rep(400, 16), n.mar = 185, n.ind = 112, type = "bc", n.pheno = 6000, latent.eff = 1.5, res.var = 1, init.seed = 92387475)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chr.len}{
%%     ~~Describe \code{chr.len} here~~
}
  \item{n.mar}{
%%     ~~Describe \code{n.mar} here~~
}
  \item{n.ind}{
%%     ~~Describe \code{n.ind} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{n.pheno}{
%%     ~~Describe \code{n.pheno} here~~
}
  \item{latent.eff}{
%%     ~~Describe \code{latent.eff} here~~
}
  \item{res.var}{
%%     ~~Describe \code{res.var} here~~
}
  \item{init.seed}{
%%     ~~Describe \code{init.seed} here~~
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
function (chr.len = rep(400, 16), n.mar = 185, n.ind = 112, type = "bc", 
    n.pheno = 6000, latent.eff = 1.5, res.var = 1, init.seed = 92387475) 
{
    set.seed(init.seed)
    mymap <- sim.map(len = chr.len, n.mar = n.mar, include.x = FALSE, 
        eq.spacing = TRUE)
    scross <- sim.cross(map = mymap, n.ind = n.ind, type = type)
    scross <- calc.genoprob(scross, step = 0)
    sim.null.pheno.data(cross = scross, n.pheno = n.pheno, latent.eff = latent.eff, 
        res.var = res.var)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
