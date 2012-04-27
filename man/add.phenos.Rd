\name{add.phenos}
\alias{add.phenos}
\title{
Add phenotypes to cross object.
}
\description{
Add phenotypes to cross object by checking index.
}
\usage{
add.phenos(cross, newdata = NULL, index = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
%%     ~~Describe \code{cross} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{index}{
%%     ~~Describe \code{index} here~~
}
}
\details{
The name \code{index} must be a phenotype in the \code{cross}
object. The row names of \code{newdata} are matched with values of \code{index}.
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
function (cross, newdata = NULL, index = NULL) 
{
    if (!is.null(newdata)) {
        if (any(names(newdata) \%in\% names(cross$pheno))) 
            warning("some cross phenotypes overwritten with new data")
        if (is.null(index)) {
            n.ind <- nind(cross)
            if (nrow(newdata) != n.ind) 
                stop(paste("newdata must have number of rows (", 
                  nrow(newdata), ") as cross individuals (", 
                  n.ind, ")", sep = ""))
            for (i in names(newdata)) cross$pheno[[i]] <- newdata[[i]]
        }
        else {
            mat <- match(as.character(cross$pheno[[index]]), 
                row.names(newdata))
            if (length(mat[is.na(mat)]) == length(mat)) 
                stop("no row names of newdata match index")
            cross$pheno <- cbind(cross$pheno, newdata[mat, , 
                drop = FALSE])
        }
    }
    cross
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
