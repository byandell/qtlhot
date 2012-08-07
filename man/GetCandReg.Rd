\name{GetCandReg}
\alias{GetCandReg}
\alias{GetCisCandReg}
\alias{GetCoMappingTraits}
\title{
Get genetic information on candidate regulators and co-mapping traits.
}
\description{
Get chromosome (phys.chr) and physical position in cM (phys.pos), along
with the LOD score (peak.lod) at the peak position (peak.pos), and the
chromosome where the peak is located (peak.chr). Some candidates may map to the
same chromosome where they are physically located.
}
\usage{
GetCandReg(cross, highlod, annot, traits, all.traits = names(cross$pheno),
  chr.pos = get.chr.pos(cross))
GetCisCandReg(cross, highlod, cand.reg, all.trait.nms = names(cross$pheno),
  chr.pos = get.chr.pos(cross))
GetCoMappingTraits(cand.reg, scan, lod.thr, drop = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cross}{
object of class \code{cross}; see \code{\link[qtl]{read.cross}}
}
  \item{highlod}{
data frame from \code{\link{pull.highlods}}, which is sparse summary of
high LODs in large \code{\link[qtl]{scanone}} object
}
  \item{cand.reg}{
data frame with candidate regulator; see value section below
}
  \item{annot}{
data frame with annotation information; must have first column as unique
identifier, third column as chromosome, and fifth column as position in
cM; typically column 2 has gene name, and column 4 has position in Mb
}
  \item{traits}{
names of traits to examine as candidate regulators; names must
correspond to phenotypes in \code{cross} object
}
  \item{lod.thr}{
critical value for LOD threshold
}
  \item{drop}{
allowably drop in LOD for support interval (default is 1.5)
}
  \item{chr.pos}{
first two columns of \code{\link[qtl]{scanone}} object with chromosome
and position in cM
}
}
\details{
Traits that map to positions close to their physical locations are said
to map in cis (local linkages).
Traits that map to positions away from their physical locations are said to map in
trans (distal linkages). There is no unambiguous way to determine how close a trait needs to
map to its physical location in order to be classified as cis. Our choice is to classify a trait as
cis if the 1.5-LOD support interval (Manichaikul et al. 2006) around the LOD peak contains
the trait’s physical location, and if the LOD score at its physical location is higher the the LOD
threshold. The function \code{GetCisCandReg} determines which of the candidate regulators map in
cis. The function \code{GetCoMappingTraits} returns a list with the putative
targets of each trait. A trait is included in the putative target list of
a trait when its LOD peak is greater than \code{lod.thr} and the
\code{drop} LOD support interval around the peak contains the location
of the trait’s QTL.
}
\value{
  \code{GetCoMappingTraits} returns a list with each element being the
  names of co-mapping traits for a particular name in \code{traits}.
  \code{GetCandReg} returns a data frame while \code{GetCisCandReg}
  returns a list with a similar candidate regulator data frame as the
  element \code{cis.reg}, and the index of trait names as the element
  \code{cis.index}. The elements of the candidate regulator data frame
  are as follows (\code{peak.pos.lower} and \code{peak.pos.upper} only
  for \code{GetCisCandReg}):
  \item{gene}{name of trait, which might be a gene name}
  \item{phys.chr}{chromosome on which gene physically resides}
  \item{phys.pos}{physical position (in cM)}
  \item{peak.chr}{chromosome where peak LOD is located}
  \item{peak.pos}{position of peak (in cM)}
  \item{peak.lod}{LOD value at peak}
  \item{peak.pos.lower,peak.pos.upper}{lower and upper bounds
    of the 1.5-LOD support interval around \code{peak.pos}}
}
\references{
Manichaikul et al. (2006) Genetics
}
\author{
Elias Chaibub Neto
}
\seealso{
  \code{\link[qtl]{scanone}},
  \code{\link{GetCoMappingTraits}}
}
\examples{
\dontrun{
## scan = scanone object with all traits
## annot = data frame with annotation of physical position
cand.reg <- GetCandReg(scan = scan, annot = annot,
  traits = names(ko.list.all), lod.thr = 3.47, drop = 1.5)
cand.reg[1:10,]
cis.cand.reg <- GetCisCandReg(cand.reg, scan, drop = 1.5, lod.thr = 3.47)
cis.cand.reg[[1]][1:10,]
comap.targets <- GetCoMappingTraits(traits = cand.reg, scan = scan,
  lod.thr = 3.47, drop = 1.5)
summary(unlist(lapply(comap.targets, length)))
}
}
\keyword{utilities}
