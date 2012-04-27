######################################################################
# perm.R
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
# Contains: NL.N.permutations, NL.N.perm, count.thr, exceed.thr, NL.N.summary
######################################################################

## This set of functions compute the permutation LOD thresholds for the NL-method
## and the permutation hotspot size thresholds for the N-method. The output
## is a list with two elements: the NL- and the N-method's threshold matrices.
## The NL-method output is a nN (number of spurious hotspot sizes) by nalpha
## (number of significance levels) threshold matrix. Note that for the
## NL-method we have a single "alpha" since we use the same significance level
## for QTL mapping and permutation significance.
## The N-method output is a nlod (number of LOD thresholds) by nalpha (number
## of significance levels) threshold matrix. Note that here we have two
## "alphas", one for the QTL mapping (the LOD thresholds) and one for the 
## permutation significance (alpha levels). 
##
NL.N.permutations <- function(cross, Ns, n.perm, alpha.levels, lod.thrs, drop=1.5,
                              verbose = FALSE)
{
  ## Time consuming step.
  mycat("NL.N.permutations ...", verbose)

  Nmax <- length(Ns)

  mats <- NL.N.perm(cross, Nmax, n.perm, lod.thrs, drop = drop,
                    verbose = verbose)
  
  NL.N.summary(mats$max.lod.quant, mats$max.N, alpha.levels)
}
NL.N.perm <- function(cross, Nmax, n.perm, lod.thrs, drop = 1.5,
                      verbose = FALSE, init.seed = 0) 
{
  set.seed(init.seed)
  n.phe <- nphe(cross)
  n.ind <- nind(cross)
  Nmax <- min(Nmax, n.phe)
  
  Ns <- seq(Nmax)
  quants <- 1 - (Ns - 1) / n.phe
  n.lod <- length(lod.thrs)

  max.N <- matrix(NA, n.perm, n.lod)
  dimnames(max.N) <- list(NULL, as.character(lod.thrs))

  max.lod.quant <- matrix(NA, n.perm, Nmax)
  dimnames(max.lod.quant) <- list(NULL, as.character(Ns))

  for(i in 1:n.perm){
    if(verbose)
      print(i)
    
    ## permute rows of the phenotype data matrix
    perm.cross <- cross
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    ## perform mapping analysis in the permuted data
    mycat("scanone", verbose)
    ## NB: scanone groups phenos in batches based on missing data patterns.
    scanmat <- scanone(perm.cross, pheno.col = c(1:n.phe), method = "hk")
    chr <- scanmat[, 1]
    scanmat <- as.matrix(scanmat[, -(1:2), drop = FALSE])

    ## apply LOD drop interval to per.scan
    mycat("LOD drop", verbose)
    scanmat <- set.to.zero.beyond.drop.int(chr, scanmat, min(lod.thrs), drop)

    ## get lod quantiles at each genomic location
    ## rows indexes the quants associated with the Ns
    ## columns indexes the genomic positions
    mycat("quantile", verbose)
    quant <- apply(scanmat, 1, quantile, quants)

    ## get the maximum lod-quantile across the genome
    ## rows indexes the permutations
    ## columns indexes the Ns
    mycat("quant", verbose)
    max.lod.quant[i,] <- apply(quant, 1, max)

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("sum", verbose)
    max.N[i,] <- apply(count.thr(scanmat, lod.thrs, FALSE), 1, max)
  }
  list(max.lod.quant = max.lod.quant, max.N = max.N)
}
count.thr <- function(scan, lod.thrs, droptwo = TRUE)
{
  ## Count number of traits at or above each value of lod.thrs for each locus.
  ## Result is n.loci * n.thr matrix.
  if(droptwo)
    scan <- scan[, -c(1,2), drop = FALSE]
  apply(scan, 1, exceed.thr, lod.thrs)
}
exceed.thr <- function(x, y)
{
  ## Finds out how many in x exceed each value of y.
  res <- rep(0, length(y))
  for(k in order(y)) {
    x <- x[x > y[k]]
    if(length(x) == 0)
      break
    res[k] <- length(x)
  }
  res
}
NL.N.summary <- function(max.lod.quant, max.N, alpha.levels)
{
  ## get the NL-method thresholds matrix
  ## rows indexes show the spurious hotspot sizes
  ## columns indexes show the QTL mapping/NL-method alphas
  Ns <- dimnames(max.lod.quant)[[2]]
  NL.thrs <- t(apply(max.lod.quant, 2, quantile, 1 - alpha.levels))
  dimnames(NL.thrs) <- list(as.character(Ns), as.character(alpha.levels))

  ## Get the N-method thresholds matrix.
  ## Entry is alpha quantile of largest hotspot size across permutations.
  ## Row indexes QTL mapping alphas; column indexes N-method alphas.
  lod.thrs <- dimnames(max.N)[[2]]
  N.thrs <- t(apply(max.N, 2, quantile, 1 - alpha.levels))
  dimnames(N.thrs) <- list(as.character(lod.thrs), as.character(alpha.levels))
  list(NL.thrs=NL.thrs, N.thrs=N.thrs)
}
