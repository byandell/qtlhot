######################################################################
# neqtl.R
#
# Elias Chaibub Neto
# Brian S Yandell
# Aimee Teo Broman
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
# Contains: pull.highlods, sexbatch.covar, scanone.permutations,
#           cat.scanone, get.tails, lod.quantile.permutation,
#           make.max.N, make.maxlod, smooth.neqtl, lod.quantile.permutations.2
######################################################################

## see  ~/p/private/diabetes1/diabetes10/scan.perm/Rcode files
pull.highlods <- function(scans, pheno.col, lod=4.5, drop.lod = 1.5,
                          extend = TRUE, restrict.lod = TRUE)
{
  if(missing(pheno.col))
    pheno.col <- names(scans)[-(1:2)]
  
  ## Extract matrix of lod scores.
  x <- as.matrix(scans[,-(1:2), drop = FALSE])

  ## Keep only traits with some LOD above lod threshold.
  keep <- apply(x, 2, function(x, lod) any(x >= lod), lod)
  x <- x[, keep, drop = FALSE]
  
  ## Find which values are at within drop.lod of maximum per chr and trait.
  if(restrict.lod) {
    ## Restrict to loci above LOD threshold.
    if(extend)
      tmpfn <- function(x, lod, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod) & (maxx < x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          (g | (c(d,0) == 1) | (c(0,d) == -1)) & (x >= lod)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod, drop.lod) {
        (max(x) <= x + drop.lod) & (x >= lod)
      }
  }
  else {
    if(extend)
      tmpfn <- function(x, lod, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod) & (maxx < x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          g | (c(d,0) == 1) | (c(0,d) == -1)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod, drop.lod) {
        maxx <- max(x)
        (maxx >= lod) & (maxx <= x + drop.lod)
      }
  }
  lodint.pos <- function(x, chr, lod, drop.lod) {
    unlist(tapply(x, chr, tmpfn, lod, drop.lod))
  }
  wh <- apply(x, 2, lodint.pos, scans$chr, lod, drop.lod)
  
  ## Get row and column indices.
  rr <- row(x)[wh]
  cc <- seq(keep)[keep][col(x)[wh]]

  ## Find which are within drop.lod of max lod per chr.
  lod <- x[wh]
  
  ## return data frame with genome row, trait column and lod value.
  out <- list(highlod = cbind.data.frame(row = rr, phenos = pheno.col[cc], lod = lod),
              chr.pos = scans[,1:2])
  class(out) <- c("highlod", "list")
  out
}

sexbatch.covar <- function(cross, batch.effect, verbose = FALSE)
{
  ic <- getsex(cross)$sex
  ## Drop sex if only one present.
  if(length(unique(ic)) == 1)
    ic <- NULL
  
  if(!is.null(batch.effect)){
    batch <- cross$pheno[,batch.effect, drop = FALSE]
    tmp <- formula(paste("~ factor(", batch.effect, ")"))
    if(verbose)
      cat("sexbatch.covar", names(tmp), levels(factor(batch[[1]])), "\n")
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    batch <- model.matrix(tmp,batch)[,-1, drop = FALSE]
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    ac <- cbind(batch,ic)
  }
  else
    ac <- ic
  
  list(addcovar = ac, intcovar = ic)
}
## Performs and saves scanone on permuted dataset
scanone.permutations <- function(cross, pheno.col = seq(3, nphe(cross)),
                                 n.perm, seed=123456789, batch.effect = NULL,
                                 pheno.set = 1,
                                 lod.min, drop.lod = 1.5)
{
  set.seed(seed[[1]])

  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  covars <- sexbatch.covar(cross, batch.effect)

  n.ind <- nind(cross)

  perms <- matrix(NA, n.ind, n.perm)

  for(i in 1:n.perm){
    perm.cross <- cross
    perms[,i] <- tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar)

    per.scan.hl <- pull.highlods(per.scan, lod = lod.min, drop.lod = drop.lod,
                                 restrict.lod = TRUE)$highlod

    save(per.scan.hl, perms,
         file=paste("per.scan",pheno.set, i,"RData",sep="."))
  }
}

## Folder should contain scanone highlods data across all traits for ONE permutation
cat.scanone <- function(dirpath = ".", filenames = permfiles, chr.pos)
{
  permfiles <- list.files(dirpath, paste("per.scan", "*", "RData", sep = "."))

  ## Make and remove per.scan.hl. Below use version from attached files.
  per.scan.hl <- NULL
  rm(per.scan.hl)
  
  for(i in 1:length(filenames)){
    attach(filenames[i], warn.conflicts = FALSE)
    if(i==1)  cat.scan.hl <- per.scan.hl else
    cat.scan.hl <- rbind.data.frame(cat.scan.hl, per.scan.hl)
    detach()
  }
  cbind.data.frame(chr.pos[cat.scan.hl$row,],cat.scan.hl)
}

## Get the 2000 highest values above lod=4 (see scan.perm.R)
get.tails <- function(highs, n.quant = 2000)
{
  n.quant <- min(n.quant, max(table(highs[,"row"])))
  tmpfn <- function(x, n) {
    l <- length(x)
    if(l < n)
      x <- c(x, rep(NA, n - l))
    rev(sort(x, na.last = FALSE))[seq(n)]
  }
  out <- tapply(highs[,"lod"], highs[,"row"], tmpfn, n.quant)
  ## Turn list of items of length n.quant into matrix.
  rows <- names(out)
  out <- matrix(unlist(out), n.quant)
  dimnames(out) <- list(seq(n.quant), rows)
  t(out)
}


## cat.scan.hl = highlods() data.frame across all traits in a tissue
## N = number of highest lod scores to save 2000 ~= 0.05 percentile
lod.quantile.permutation <- function(cat.scan.hl,N,lod.thr,window,chr.pos,n.phe)
{
  n.chr <- levels(chr.pos$chr)

  ## Elias quantiles.
  quant <- get.tails(cat.scan.hl, n.quant = N)
  N <- ncol(quant)
  max.lod.quant <- apply(quant,2,max,na.rm=TRUE)
  names(max.lod.quant) <- paste(paste(round(1-as.numeric(dimnames(quant)[[2]])/n.phe,4)*100,
                                      "%", sep=""),1:N, sep="_")

  ## Jansen counts.
  ## Count number of LODs per position and find max (raw or smoothed).
  max.N <- make.max.N(cat.scan.hl, lod.thr, chr.pos, window)
  
  list(max.lod.quant=max.lod.quant, max.N=max.N)
}
make.max.N <- function(cat.scan.hl, lod.thr, chr.pos, window,
                       max.hl = make.maxlod(cat.scan.hl, chr.pos))
{
  l.lod.thr <- length(lod.thr)
  max.N <- data.frame(max.N = vector(length = l.lod.thr), max.N.win = vector(length = l.lod.thr),
                      row.names = lod.thr, check.names = TRUE)
  for(j in 1:l.lod.thr){
    XX <- cat.scan.hl$lod >= lod.thr[j]
    max.N$max.N[j] <- max(tapply(XX, cat.scan.hl$row, sum, na.rm = TRUE), na.rm = TRUE)

    neqtl.pos <- smooth.neqtl(cat.scan.hl, chr.pos, max.hl, lod.thr[j], window)
    
    max.N$max.N.win[j] <- max(neqtl.pos[,3])
  }
  max.N
}
make.maxlod <- function(cat.scan.hl, chr.pos)
{
  ## find high LOD and position per chromosome.
  n.chr <- levels(chr.pos$chr)

  tmpfn <- function(x) {
    if(is.null(x))
      0
    else
      max(x, na.rm = TRUE)
  }
  tmpfn2 <- function(x) {
    if(is.null(x))
      NA
    else
      mean(x, na.rm=TRUE)
  }
  tmpfn3 <- function(a) {
    is.nan(a) | a==max(a, na.rm=TRUE)
  }
  
  maxlod.hl <- maxlod.pos.hl <- vector("list", length(n.chr))
  names(maxlod.pos.hl) <- n.chr
  for(k in seq(along=n.chr)) {
    scan.out.bychr <- cat.scan.hl[cat.scan.hl$chr==n.chr[k],]
    ## This is kludgey. How to make more efficient?
    scan.out.bychr$phenos <- ordered(scan.out.bychr$phenos, unique(scan.out.bychr$phenos))
    tmp <- tapply(scan.out.bychr$lod, scan.out.bychr$phenos, tmpfn3)
    scan.out.bychr <- scan.out.bychr[unlist(tmp),]
    if(nrow(scan.out.bychr)) {
      scan.out.bychr$phenos <- ordered(scan.out.bychr$phenos, unique(scan.out.bychr$phenos))
      ## Find high lod.
      maxlod.hl[[k]] <- tapply(scan.out.bychr$lod, scan.out.bychr$phenos, tmpfn)
      ## Find position of high lod.
      maxlod.pos.hl[[k]] <- tapply(scan.out.bychr$pos, scan.out.bychr$phenos, tmpfn2)
    }
  }
  list(lod = maxlod.hl, pos = maxlod.pos.hl)
}

smooth.neqtl <- function(cat.scan.hl, chr.pos, max.hl = make.maxlod(cat.scan.hl, chr.pos),
                         lod.thr, window = 5)
{
  chr <- chr.pos$chr
  pos <- chr.pos$pos
  n.chr <- levels(chr.pos$chr)

  maxlod.thr.pos <- max.hl$pos
  for(k in seq(along=n.chr))
    maxlod.thr.pos[[k]] <- max.hl$pos[[k]][max.hl$lod[[k]] >= lod.thr]
  
  out <- smoothall(maxlod.thr.pos,thechr = chr.pos$chr, thepos = chr.pos$pos, window = window)
  ## Recover marker information.
  rownames(out) <- rownames(chr.pos)
  out
}

## Run permuations
lod.quantile.permutations.2 <- function(cross, pheno.col, N, n.perm, lod.thr, 
                                        batch.effect, window, seed = 123456789, verbose = FALSE)
{

  ## Set up pseudo-random number generation seeds.
  set.seed(seed[[1]])
  all.seeds <- sample(c(98765:987654), n.perm, replace=FALSE)

  ## Set up matrices to record values.
  nphe.cross <- nphe(cross)
  N <- N[N <= nphe.cross]
  quants <- 1 - (N - 1)/nphe.cross
  l.N <- length(N)
  max.lod.quant <- matrix(NA, n.perm, l.N)
  dimnames(max.lod.quant)[[2]] <- paste(paste(round(quants,4)*100, "%", sep=""), 
                                        N, sep="_")
  l.lod.thr <- length(lod.thr)
  max.N <- matrix(NA, n.perm, l.lod.thr)
  max.lod.quant <- matrix(NA, n.perm, l.N)
  max.N.window <- matrix(NA, n.perm, l.lod.thr)

  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  covars <- sexbatch.covar(cross, batch.effect)

  n.ind <- nind(cross)

  for(i in 1:n.perm){
    perm.cross <- cross
    set.seed(all.seeds[i])
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar)

    ## Elias' quantiles.
    quant <- apply(per.scan[,-(1:2)], 1, quantile, quants)
    max.lod.quant[i,] <- apply(quant,1,max)

    ## Jansen's count.
    max.N[i,] <- apply(count.thr(per.scan, lod.thr, droptwo = TRUE), 1, sum)

    ## Smoothed count.
    maxlod <-  apply(per.scan[,-(1:2)], 2, tapply, per.scan[,1], max)
    chrs <- dimnames(maxlod)[[1]]
    chr <- factor(per.scan[,1], levels = chrs)
    pos <- as.numeric(per.scan[,2])

    maxlod.pos <- maxlod.thr.pos <- vector("list", length(chrs))
    names(maxlod.pos) <- names(maxlod.thr.pos) <- chrs
    for(k in seq(length(chrs))) {
      scan.out.bychr <- per.scan[per.scan[,1] == chrs[k], ]
      maxlod.pos[[k]] <- apply(scan.out.bychr[,-(1:2)], 2, function(a,b)
                               mean(b[!is.nan(a) & a==max(a, na.rm=TRUE)]), scan.out.bychr[,2])
    }

    for(j in seq(along = lod.thr)){
      for(k in chrs)
        maxlod.thr.pos[[k]] <- maxlod.pos[[k]][maxlod[k,] >= lod.thr[j]]
      neqtl.pos <- smoothall(maxlod.thr.pos, thechr = chr, thepos = pos, window = window)
      max.N.window[i,j] <- max(neqtl.pos[,3])
    }
    print(i)
  }
  list(max.lod.quant=max.lod.quant, 
       max.N=max.N,
       max.N.window=max.N.window)
}
