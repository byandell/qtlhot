\name{PrecTpFpMatrix}
\alias{PrecTpFpMatrix}
\title{
Determine false positive and true positive rates for known targets.
}
\description{
Determine how well different tests do to predict candidates of regulation.
}
\usage{
PrecTpFpMatrix(alpha, nms, val.targets, all.orfs, to.load, cis.index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
significance levels at which summaries are computed
}
  \item{nms}{
names of traits that are candidate regulators
}
  \item{val.targets}{
validated targets of candidate regulators
}
  \item{all.orfs}{
all trait names
}
  \item{to.load}{
name of file with RData object of joined output
}
  \item{cis.index}{
numeric index into \code{all.orfs} to identify cis traits
}
}
\details{
After loading the joined results we use the function PrecTpFpMatrix to summarize the
performance of the different methods in terms of ``biologically validated'' true positives, false
positives and precision, of the inferred causal relations. We define a
true positive as a statistically significant causal relation between a
gene and a putative target gene when the putative target gene belongs to
the known signature of the gene. Similarly, we define a false positive
as a statistically significant causal relation between a gene and a
putative target gene when the target gene doesn’t belong to the
signature. (For the AIC and BIC methods that do not provide a p-value
measuring the significance of the causal call, we simply use the
detected causal relations in the computation of true and false
positives). The “validated precision”, is computed as the ratio of true
positives by the sum of true and false positives. The PrecTpFpMatrix
computes these measures to both all ko-genes, and to cis ko-genes
only.
}
\value{
  List containing
  \item{Prec1,Prec2}{matrix of precision with rows for significance level and
    columns for test; first is for all, second is for cis candidates only}
  \item{Tp1,Tp2}{matrix of true positive rate with rows for significance level and
    columns for test; first is for all, second is for cis candidates only}
  \item{Fp1,Fp2}{matrix of false positive rate with rows for significance level and
    columns for test; first is for all, second is for cis candidates only}
}
\author{
Elias Chaibub Neto
}
\seealso{
\code{\link{GetCandReg}}
}
\examples{
\dontrun{
aux <- PrecTpFpMatrix(alpha = seq(0.01, 0.10, by = 0.01),
  nms = cand.reg[, 1], val.targets = ko.list,
  all.orfs = names(yeast.orf$pheno),
  to.load = "joined_ko_output.RData",
  cis.index = cis.cand.reg$cis.index)
}
}
\keyword{utilities}
