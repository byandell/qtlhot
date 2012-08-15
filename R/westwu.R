######################################################################
# westwu.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: WW.permutations, WW.perm, WW.summary
######################################################################

## This function computes the West/Wu permutation thresholds.
## The output is a nlod (number of LOD thresholds) by nalpha (number of 
## significance levels) matrix, where each entry shows the hotspot size 
## significance threshold of the West/Wu approach. Note we have two "alphas" 
## here, one for the QTL mapping (the LOD thresholds) and one for the 
## permutation significance (alpha levels).
##
## Note that I separated the original WW.permutations() into a piece that do the 
## actual permutations [WW.perm.matrix() function] and a piece that summarizes it
## [the WW.summary() function] in the same way you did with the NL.N.permutations()
## function.  
## 
WW.permutations <- function(scanmat, lod.thrs, alpha.levels, n.perm, verbose = FALSE) 
{
  ## Time consuming step.
  mycat("WW.permutations ...", verbose)

  max.WW <- WW.perm(scanmat, lod.thrs, n.perm, verbose)
  WW.summary(max.WW, alpha.levels)
}

WW.perm <- function(scanmat, lod.thrs, n.perm, verbose = FALSE)
{
  ## Permute all but first column.
  nlod <- length(lod.thrs)
  max.WW <- matrix(NA, n.perm, nlod)
  dimnames(max.WW) <- list(NULL, as.character(lod.thrs))
  
  for(i in 1:n.perm) {
    if(verbose)
      print(i)

    ## Separately permute columns of the scan object separately by trait.
    mycat("sample", verbose)
    scanmat <- apply(scanmat, 2, sample)

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("sum", verbose)
    max.WW[i,] <- apply(count.thr(scanmat, lod.thrs, FALSE), 1, max)
  }
  max.WW
}

WW.perm.highlod <- function(highobj, lod.thrs, n.perm, verbose = FALSE)
{
  phenos <- sort(unique(highobj$highlod$pheno))
  n.chrpos <- nrow(highobj$chr.pos)

  ## Want to replace row in highobj$highlod with permuted row.
  ## With separate permutation by pheno of seq(n.chrpos).
  
}

WW.summary <- function(max.WW, alpha.levels)
{
  nalpha <- length(alpha.levels)
  WW.thrs <- t(apply(max.WW, 2, quantile, 1 - alpha.levels))
  dimnames(WW.thrs) <- list(dimnames(max.WW)[[2]], as.character(alpha.levels))
  WW.thrs
}




