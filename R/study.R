######################################################################
# study.R
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
# Contains: mySimulations, get.hotspot, sim.hotspot, filter.threshold,
#           NL.counts, N.WW.counts, mycat
######################################################################

## This function computes the error rates for the NL-, N- and West/Wu methods
## out of nSim simulations. (It also returns the method's thresholds).
## At each simulation iteration this function: (1) generates a null dataset 
## cross; (2) perform Haley-Knott mapping for all traits; (3) applies the LOD
## drop interval computation to the scanone object; (4) determine the NL-, N-
## and West/Wu approches thresholds using the scanone results from step 3; 
## (5) for each of the three methods it computes the proportion of times out 
## of the nSim simulations we detected at least one false hotspot anywhere in
## the genome.  
##
mySimulations <- function(...) sim.hotspot(...)
get.hotspot <- function(filenames,
                        ## Following supplied in filenames[1].
                        Nmax, out.sim)
{
  ## filenames = list.files(".", paste(prefix, latent.eff, sets, "RData", sep = "."))
  ## latent.eff = 0, prefix = "Pilot", sets = "[0-9][0-9]*"
  
  ## Get stored Nmax value, and out.sim to get alpha.levels and lod.thrs.
  load(filenames[1])
  nSim <- length(filenames)
  Nmax <- Nmax ## Null action to make sure Nmax is valued.
  Ns <- seq(Nmax)
  tmp <- dimnames(out.sim$N.thrs)
  lod.thrs <- as.numeric(tmp[[1]])
  alpha.levels <- as.numeric(tmp[[2]])

  ## May not have names on rows and columns.
  nalpha <- ncol(out.sim$N.thrs)
  nlod <- nrow(out.sim$N.thrs)
  
  ## outputs count the number of times we detected
  ## a hotspot using the respective method
  outNL <- matrix(0, Nmax, nalpha)
  dimnames(outNL) <- list(NULL, alpha.levels)
  outN <- outWW <- matrix(0, nlod, nalpha)
  dimnames(outN) <- dimnames(outWW) <- list(lod.thrs, alpha.levels)
  
  ## we are saving the thresholds of each simulation
  thrNL <- array(dim=c(Nmax, nalpha, nSim))
  thrN <- thrWW <- array(dim=c(nlod, nalpha, nSim))
  dimnames(thrN) <- dimnames(thrWW) <- list(lod.thrs, alpha.levels, NULL)
  dimnames(thrNL) <- list(NULL, alpha.levels, NULL)
  
  for(k in 1:nSim) {
    mycat(k, TRUE, TRUE)
    load(filenames[k])

    thrNL[,,k] <- out.sim$NL.thrs
    thrN[,,k] <- out.sim$N.thrs
    thrWW[,,k] <- out.sim$WW.thrs    
    outNL <- outNL + out.sim$NL
    outN <- outN + out.sim$N.counts
    outWW <- outWW + out.sim$WW.counts
  }
  
  NL.err <- outNL/nSim
  dimnames(NL.err) <- list(as.factor(Ns), as.factor(alpha.levels))
  N.err <- outN / nSim
  dimnames(N.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  WW.err <- outWW / nSim
  dimnames(WW.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  list(nSim = nSim, NL.err=NL.err, N.err=N.err, WW.err=WW.err, thrNL=thrNL, thrN=thrN, 
       thrWW=thrWW)  
}
sim.hotspot <- function(nSim, 
                        cross, 
                        nT,
                        latent.eff,
                        res.var = 1,
                        Ns,
                        n.perm,
                        alpha.levels,
                        lod.thrs,
                        drop=1.5,
                        verbose = FALSE)
{
  Nmax <- length(Ns)

  nalpha <- length(alpha.levels)
  nlod <- length(lod.thrs)

  ## outputs count the number of times we detected
  ## a hotspot using the respective method
  outNL <- matrix(0, Nmax, nalpha)
  outN <- outWW <- matrix(0, nlod, nalpha)

  ## we are saving the thresholds of each simulation
  thrNL <- array(dim=c(Nmax, nalpha, nSim))
  thrN <- array(dim=c(nlod, nalpha, nSim))
  thrWW <- array(dim=c(nlod, nalpha, nSim))

  for(k in 1:nSim){
    mycat(k, verbose, TRUE)

    mycat("sim.null.pheno.data", verbose)
    ncross <- sim.null.pheno.data(cross, nT, latent.eff, res.var)
  
    ## Simulate correlated phenotypes and create threshold summaries.
    out.sim <- filter.threshold(ncross, nT, latent.eff[k], res.var,
                             lod.thrs, drop,
                             Ns, n.perm, alpha.levels,
                             verbose)

    thrNL[,,k] <- out.sim$NL.thrs
    thrN[,,k] <- out.sim$N.thrs
    thrWW[,,k] <- out.sim$WW.thrs    
    outNL <- outNL + out.sim$NL
    outN <- outN + out.sim$N.counts
    outWW <- outWW + out.sim$WW.counts
  }

  
  NL.err <- outNL/nSim
  dimnames(NL.err) <- list(as.factor(Ns), as.factor(alpha.levels))
  N.err <- outN / nSim
  dimnames(N.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  WW.err <- outWW / nSim
  dimnames(WW.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  list(nSim = nSim, NL.err=NL.err, N.err=N.err, WW.err=WW.err, thrNL=thrNL, thrN=thrN, 
       thrWW=thrWW)  
}
########################################################################################
filter.threshold <- function(cross, nT, latent.eff, res.var,
                             lod.thrs, drop = 1.5,
                             Ns, n.perm, alpha.levels,
                             NL.N.thrs = NL.N.permutations(cross, Ns, n.perm, alpha.levels,
                               lod.thrs, verbose = verbose),
                             WW.thrs = WW.permutations(scan.drop, lod.thrs, alpha.levels, n.perm),
                             verbose = FALSE)
{
  mycat("scanone", verbose)
  scanmat <- scanone(cross, pheno.col=c(1:nT), method="hk")
  chr <- scanmat[,1]
  scanmat <- as.matrix(scanmat[, -(1:2), drop = FALSE])
  
  ## we use the lowest lod threshold in the LOD drop interval 
  mycat("set.to.zero.beyond.drop.int", verbose)
  scan.drop <- set.to.zero.beyond.drop.int(chr, scanmat, min(lod.thrs), drop)

  ## computes an array of size Nmax by nalpha by npos.
  ## showing for each Ns size and alpha level, the 
  ## hotspot sizes at each genomic location.
  mycat("NL.counts", verbose)
  NL.thrs <- NL.N.thrs[[1]]
  N.thrs <- NL.N.thrs[[2]]
  Nmax <- length(Ns)
  NL <- NL.counts(scan.drop, Nmax, NL.thrs)
  
  ## computes a matrix of size nlod by npos.
  ## showing for each lod threshold the 
  ## hotspot sizes at each genomic location.
  mycat("N.WW.counts", verbose)
  N.WW <- N.WW.counts(scan.drop, lod.thrs, N.thrs, WW.thrs)
  
  list(NL.thrs = NL.thrs, N.thrs = N.thrs, WW.thrs = WW.thrs, NL = NL,
       N.counts = N.WW$N, WW.counts = N.WW$WW)
}

## Computes an array of size Nmax (number of spurious hotspots sizes) by 
## nalpha (number of significance levels) by npos (number of locus), and for
## each spurious hotspot size/significance level threshold, it computes the 
## number of traits mapping with LOD higher than the threshold at each one
## of the genomic positions.
##
NL.counts <- function(scanmat, Nmax, NL.thrs)
{
  ## get the maximum spurious hotspot size (N-method) 
  ## for different QTL mapping significance levels

  n.phe <- ncol(scanmat)

  quants <- 1 - (seq(Nmax) - 1) / n.phe
  XX <- apply(apply(scanmat, 1, quantile, quants), 1, max)
  NL.counts <- apply(NL.thrs, 2, function(x,y) (x < y), XX)
  ## dimnames(NL.counts)[[2]] <- seq(Nmax)
  NL.counts
}

## Computes a matrix of size nlod (number of mapping thresholds) by npos 
## (number of locus), and for each LOD threshold, it computes the number
## of traits mapping with LOD higher than the threshold at each one of
## the genomic positions. The same counts are used by the N- and WW-methods.
##
N.WW.counts <- function(scanmat, lod.thrs, N.thrs, WW.thrs)
{
  ## XX = genome position by number of traits above LOD threshold.
  XX <- apply(count.thr(scanmat, lod.thrs, FALSE), 1, max)

  ## N.counts[lod,alpha] = TRUE if max hotspot size using lod is above alpha perm threshold.
  N.counts <- apply(N.thrs, 2, function(x,y) (x < y), XX)

  ## WW.counts[lod,alpha] = TRUE if max hotspot size using lod is above alpha perm threshold.
  WW.counts <- apply(WW.thrs, 2, function(x,y) (x < y), XX)
  dimnames(N.counts) <- dimnames(WW.counts) <- dimnames(N.thrs)

  list(N = N.counts, WW = WW.counts)
} 
mycat <- function(title, verbose = FALSE, init = FALSE)
{
  if(verbose) {
    if(verbose > 1) {
      if(init)
        cat("user system elapsed time\n")
      else
        cat(round(as.numeric(proc.time()[1:3])), "")
    }
    cat(title, "\n")
  }
}
