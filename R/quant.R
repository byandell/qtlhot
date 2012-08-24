######################################################################
# quant.R
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
# Contains: quant.lod, plot.quant.lod
# legacy: qtlhot.scan, quant.plot
######################################################################
## This file concerns the collection of permutations on multiple traits.
## quant.scanone is dubious, but quant.sum is important.

quant.scanone <- function(scan, max.lod.quant, lod.thrs,
                       probs = seq(length(lod.thrs)) * 0.01, level = 0.05)
{
  thr.level <- min(which(probs >= level))
  lod.thr <- lod.thrs[thr.level]
  
  scan.hl <- highlod(scan, lod = lod.thr, restrict.lod = TRUE)
  quant <- quant.slide(max.lod.quant, lod.thrs, probs, level = level)
  
  hotsize(scan.hl, lod.thr, quant.level = quant)
}
###################################################################################################
quant.lod <- function(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)
{
  ## lod.thrs and probs should be as provided to create summaries.
  ## Should in future embed these.
  
  thr.level <- min(which(probs >= level))
  lod.thr <- lod.thrs[thr.level]
  
  slider <- quant.slide(max.lod.quant, lod.thrs, probs, level, TRUE)

  quant.N.window <- apply(max.N.window, 2, quantile, probs = probs, na.rm = TRUE)
  quant.N <- apply(max.N, 2, quantile, probs = probs, na.rm = TRUE)
  dimnames(quant.N)[[2]] <- dimnames(quant.N.window)[[2]] <- rev(paste(lod.thr, lod.thrs))

  out <- list(quant.N.window = quant.N.window, quant.N = quant.N,
              quant = slider$quant, max.quant = slider$max.quant, quant.level = slider$quant.level,
              thr.level = thr.level, lod.thr = lod.thr)
  class(out) <- c("quant.lod", "list")
  attr(out, "lod.thrs") <- lod.thrs
  attr(out, "probs") <- probs

  out
}
###################################################################################################
plot.quant.lod <- function(x, ...)
{
  thr.level <- x$thr.level
  lod.thr <- x$lod.thr
  n.probs <- length(probs)

  tmp.plot <- function(x.vals, quant, x.crit, probs, level, is.quantile = FALSE, main = "",
                       add.level = FALSE)
  {
    n.probs <- length(probs)
    thr.level <- min(which(probs >= level))
    quant.thr <- rev(quant[thr.level,])[thr.level]

    xlabs <- "single trait LOD threshold"
    if(is.quantile)
      xlabs <- paste(xlabs, "quantile")
    
    plot(range(x.vals), c(0, max(quant)), type = "n", xlab = "", ylab = "")
    mtext(xlabs, 1, 2)
    mtext("hotspot size", 2, 2)
    abline(v = x.crit, col = "darkgray", lty = 2)
    abline(h = quant.thr, col = "darkgray", lty = 2)
    mtext(ceiling(quant.thr), 2, at = quant.thr, las = 2, cex = 0.5)
    for(i in seq(along = probs)) {
      lines(rev(sort(x.vals)), quant[i,],
            lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }
    text(x.crit, rev(quant[n.probs,])[thr.level] + 5, 1 - max(probs), adj = 0)
    text(x.crit, rev(quant[1,])[thr.level] - 5, 1 - min(probs), adj = 1)
    text(x.vals[n.probs], quant.thr + 5, 1 - level, adj = 1)
    if(add.level)
      main <- paste(main, "\n hotspot size significance level =", 1 - max(probs), "to", 1 - min(probs))
    mtext(main, 3, 0.5)

    quant.thr
  }

  tmpar <- par(mfrow = c(2,2), mar = c(3.1,3.1,2.1,0.1))
  ## Jansen method, smoothing.
  x$quant.thr <- tmp.plot(lod.thrs, x$quant.N.window, lod.thr, probs, level, FALSE,
                            "Jansen method 5cM window")
  
  tmp.plot(probs, x$quant.N.window, level, probs, level, TRUE,
           "Jansen method 5cM window")

  tmp.plot(lod.thrs, x$quant.N, lod.thr, probs, level, FALSE,
           "Jansen method per locus")
  tmp.plot(probs, x$quant.N, level, probs, level, TRUE,
           "Jansen method per locus")
  par(tmpar)
  
  ## Chaibub Neto method.
  plot(c(1,x$max.quant), c(min(x$quant, na.rm = TRUE), max(x$quant, na.rm = TRUE)), type = "n",
       xlab = "significant hotspot size with given threshold",
            ylab = "hotspot LOD score threshold",
            log = "xy")
  abline(h = lod.thr, col = "darkgray", lty = 2)
  abline(v = x$quant.thr, col = "darkgray", lty = 2)
  mtext(ceiling(x$quant.thr), 1, at = x$quant.thr, las = 2)
  for(i in seq(along = probs)) {
    tmp <- (x$quant[i,seq(x$max.quant)] >= lod.thr)
    if(any(tmp))
      lines(seq(x$max.quant)[tmp], x$quant[i,seq(x$max.quant)][tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0))
    if(any(!tmp))
      lines(seq(x$max.quant)[!tmp], x$quant[i,seq(x$max.quant)][!tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0),
            col = "darkgray")
  }
  
  n.thr2 <- length(lod.thrs) / 2
  text(n.thr2 + 1, x$quant[n.probs, n.thr2], 1 - max(probs), adj = 0)
  text(n.thr2 - 1, x$quant[1,n.thr2], 1 - min(probs), adj = 1)
  title(paste("hotspot LOD threshold by hotspot size\nsignificance level =",
              1 - max(probs), "to", 1 - min(probs)))
  
  invisible(x)
}
###################################################################################################
quant.slide <- function(max.lod.quant, lod.thrs, probs, level = 0.95, show.max = FALSE)
{
  ## level and probs must be > 0.5. If not, use 1 - x. Probably want to change this later.
  if(level > 0.5)
    level <- 1 - level
  if(min(probs) > 0.5)
    probs <- 1 - probs

  thr.level <- min(which(probs <= level))
  lod.thr <- lod.thrs[thr.level]
  
  quant <- apply(max.lod.quant, 2, quantile, probs = probs, na.rm = TRUE)
  max.quant <- sum(apply(quant, 2, function(x) !all(is.na(x))))
  quant.level <- quant[which((round(probs - level, 2) == 0)), seq(max.quant)]
  quant.level <- quant.level[quant.level >= lod.thr]

  lq <- length(quant.level)
  if(lq) {
    if(lq < max.quant & quant.level[lq] > lod.thr)
      quant.level[as.character(lq + 1)] <- lod.thr
  }
  
  if(show.max)
    list(quant = quant, max.quant = max.quant, quant.level = quant.level)
  else
    quant.level
}
