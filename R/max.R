######################################################################
# max.R
#
# Elias Chaibub Neto
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
# Contains: max.hotsize, max.highlod, quantile.highlod, get.tails
######################################################################

max.hotsize <- function(x, ...)
{
  if(is.null(x))
    return(NULL)
  
  class(x) <- class(x)[-1]
  tmpmax <- function(x, lc) max(x, lodcolumn = lc)[,c(1,2,2+lc)]
  lc <- 1
  out <- tmpmax(x, lc)
  ## max.N.window
  if(!is.null(attr(x, "window"))) {
    lc <- lc + 1
    out <- cbind(out, tmpmax(x, lc))
  }
  if(!is.null(attr(x, "quant.level"))) {
    lc <- lc + 1
    out <- cbind(out, tmpmax(x, lc))
  }
  out
}
max.highlod <- function(x, lod.thr = NULL, ...)
{
  if(length(lod.thr) > 1) {
    out <- NULL
    out.thr <- NULL
    for(lod in lod.thr) {
      tmp <- hotsize(x, lod, ...)
      if(!is.null(tmp)) {
        out <- rbind(out, max(tmp))
        out.thr <- c(out.thr, lod)
      }
    }
    if(!is.null(out))
      out$lod.thr <- out.thr
    
    class(out) <- c("summary.scanone", "data.frame")
    out
  }
  else
    max(hotsize(x, lod.thr, ...))
}
quantile.highlod <- function(x, probs = NULL, n.quant = 2000, n.phe, ...)
{
  highlod <- highlod.thr(x, lod.thr)
  lod.thr <- highlod$lod.thr
  highlod <- highlod$highlod
  if(!nrow(highlod))
    return(NULL)
  
  ## Probabilities if requested.
  if(!is.null(probs)) {
    if(missing(n.phe))
      stop("need to supply n.phe along with probs")
    s.quant <- ceiling(probs * n.phe)
    n.quant <- max(s.quant)
  }
  else
    s.quant <- seq(n.quant)
  
  ## Quantiles from 1 to n.quant.
  if(n.quant)
    get.tails(highlod, n.quant, s.quant)
  else
    NULL
}
#####################################################################
## Get the n.quant highest values above lod=4 (see scan.perm.R)
get.tails <- function(highs, n.quant = 2000, s.quant = seq(n.quant))
{
  ## Limit n.quant to range of data.
  n.quant <- min(n.quant, max(s.quant), max(table(highs[,"row"])))
  s.quant <- s.quant[s.quant <= n.quant]
  
  tmpfn <- function(x, sn) {
    x <- sort(x, decreasing = TRUE)
    x[sn]
  }
  
  out <- tapply(highs[,"lod"], highs[,"row"], tmpfn, s.quant)
  ## Turn list of items of length n.quant into matrix. This step takes time!
  out <- matrix(unlist(out), n.quant)
  out <- apply(out, 1, max, na.rm = TRUE)
  names(out) <- s.quant

  out
}
#####################################################################
quantile.hotperm <- function(x, probs = 0.95, ..., lod.thr = NULL)
{
  if(max(probs) <= 0.5)
    probs <- 1 - probs
  
  myquant <- function(x, probs) {
    x[is.na(x)] <- 0
    out <- as.matrix(apply(x, 2, quantile, probs = probs))
    if(length(probs) > 1)
      out <- t(out)
    out
  }
  
  out <- list()
  prob.names <- as.character(signif(probs, 3))

  ## max.N
  out$max.N <- myquant(x$max.N, probs)
  dimnames(out$max.N) <- list(signif(attr(x, "lod.thrs"), 3), prob.names)

  ## max.N.window
  if(!is.null(x$max.N.window)) {
    out$max.N.window <- myquant(x$max.N.window, probs)
    dimnames(out$max.N.Window) <- dimnames(out$max.N)
  }
  
  ## max.lod.quant
  if(!is.null(x$max.lod.quant)) {
    quant <- myquant(x$max.lod.quant, probs)
    if(!is.null(lod.thr)) {
      tmp <- quant <= lod.thr & quant > 0
      if(any(tmp))
        quant[tmp] <- 0
      first.zero <- apply(quant, 2,
                          function(x) {
                            z <- x == 0
                            if(any(z))
                              min(which(z))
                            else
                              0
                          })
      offset <- nrow(quant) * (seq(ncol(quant)) - 1)
      first.zero <- first.zero[first.zero>0] + offset[first.zero>0]
      quant[first.zero] <- lod.thr
    }
    dimnames(quant) <- list(as.character(dimnames(x$max.lod.quant)[[2]]), prob.names)
    quant <- quant[apply(quant, 1, function(x) any(x > 0)),, drop = FALSE]
    quant[quant == 0] <- NA
    out$max.lod.quant <- quant
  }
  out
}
