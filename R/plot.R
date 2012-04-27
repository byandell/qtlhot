######################################################################
# plot.R
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
# Contains: pull.hotspots, quant.sum, quant.plot, hotspot.scan, scan.hl.plot, hotspot.plot
######################################################################

## Routines to show hotspots. None here yet.
pull.hotspots <- function(cross, scan.hl, chr.pos = NULL, lod.thr = 5, slide.thr = NULL, verbose = FALSE)
{
  if(is.null(chr.pos)) {
    n.traits <- length(all.traits)
    cross$pheno$trait <- rnorm(nind(cross))
    chr.pos <- scanone(cross, pheno.col= find.pheno(cross, "trait"))[,1:2]
  }
  
  scan.hl <- scan.hl[scan.hl$lod >= lod.thr,]
  if(nrow(scan.hl)) {
    ## Straight count of LODs above threshold. Includes shoulders and peaks.
    tbl <- table(scan.hl$row)
    tmp <- rep(0, nrow(chr.pos))
    if(length(tbl))
      tmp[as.numeric(names(tbl))] <- tbl
    chr.pos$max.N <- tmp

    if(!is.null(slide.thr)) {
      ## If sliding thresholds supplied.
      slide.thr <- slide.thr[slide.thr >= lod.thr & !is.na(slide.thr)]
      if(length(slide.thr)) {
        ## Want to work down slide.thr. One way is to see how many are above smallest
        ## and stop if above the index.
        index <- rep(TRUE, nrow(scan.hl))
        hot.size <- rep(0, nrow(chr.pos))
        for(hot.crit in rev(seq(length(slide.thr)))) {
          if(verbose)
            cat(hot.crit, "\n")
          above <- scan.hl$lod[index] >= slide.thr[hot.crit]
          if(!any(above))
            break;
          tbl <- table(scan.hl$row[index][above])
          tbl <- tbl[tbl >= hot.crit]
          if(length(tbl)) {
            ## Record hotspot size if above hot.crit for thr, then mask those loci out.
            hot.size[as.numeric(names(tbl))] <- tbl
            index[scan.hl$row %in% names(tbl)] <- FALSE
          }
          if(!any(index))
            break;
        }
        chr.pos$quant <- hot.size
      }
    }
    chr.pos
  }
  else
    NULL
}
###################################################################################################
quant.sum <- function(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)
{
  ## lod.thrs and probs should be as provided to create summaries.
  ## Should in future embed these.
  
  thr.level <- min(which(probs >= level))
  lod.thr <- lod.thrs[thr.level]
  
  quant <- apply(max.lod.quant, 2, quantile, probs = probs, na.rm = TRUE)
  max.quant <- sum(apply(quant, 2, function(x) !all(is.na(x))))
  quant.level <- quant[which((round(probs - level, 2) == 0)), seq(max.quant)]
  quant.level <- quant.level[quant.level >= lod.thr]

  quant.N.window <- apply(max.N.window, 2, quantile, probs = probs, na.rm = TRUE)
  quant.N <- apply(max.N, 2, quantile, probs = probs, na.rm = TRUE)
  dimnames(quant.N)[[2]] <- dimnames(quant.N.window)[[2]] <- rev(paste(lod.thr, lod.thrs))

  list(quant.N.window = quant.N.window, quant.N = quant.N,
       quant = quant, max.quant = max.quant, quant.level = quant.level,
       thr.level = thr.level, lod.thr = lod.thr)
}
quant.plot <- function(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level = 0.95)
{
  ## lod.thrs and probs should be as provided to create summaries.
  ## Should in future embed these.

  out <- quant.sum(max.lod.quant, max.N, max.N.window, lod.thrs, probs, level)
  
  thr.level <- out$thr.level
  lod.thr <- out$lod.thr
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
  out$quant.thr <- tmp.plot(lod.thrs, out$quant.N.window, lod.thr, probs, level, FALSE,
                            "Jansen method 5cM window")
  
  tmp.plot(probs, out$quant.N.window, level, probs, level, TRUE,
           "Jansen method 5cM window")

  tmp.plot(lod.thrs, out$quant.N, lod.thr, probs, level, FALSE,
           "Jansen method per locus")
  tmp.plot(probs, out$quant.N, level, probs, level, TRUE,
           "Jansen method per locus")
  par(tmpar)
  
  ## Chaibub Neto method.
  plot(c(1,out$max.quant), c(min(out$quant, na.rm = TRUE), max(out$quant, na.rm = TRUE)), type = "n",
       xlab = "significant hotspot size with given threshold",
            ylab = "hotspot LOD score threshold",
            log = "xy")
  abline(h = lod.thr, col = "darkgray", lty = 2)
  abline(v = out$quant.thr, col = "darkgray", lty = 2)
  mtext(ceiling(out$quant.thr), 1, at = out$quant.thr, las = 2)
  for(i in seq(along = probs)) {
    tmp <- (out$quant[i,seq(out$max.quant)] >= lod.thr)
    if(any(tmp))
      lines(seq(out$max.quant)[tmp], out$quant[i,seq(out$max.quant)][tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0))
    if(any(!tmp))
      lines(seq(out$max.quant)[!tmp], out$quant[i,seq(out$max.quant)][!tmp],
            lwd = 1 + 2 * (round(probs[i] - 0.95, 2) == 0),
            col = "darkgray")
  }
  
  n.thr2 <- length(lod.thrs) / 2
  text(n.thr2 + 1, out$quant[n.probs, n.thr2], 1 - max(probs), adj = 0)
  text(n.thr2 - 1, out$quant[1,n.thr2], 1 - min(probs), adj = 1)
  title(paste("hotspot LOD threshold by hotspot size\nsignificance level =",
              1 - max(probs), "to", 1 - min(probs)))
  
  invisible(out)
}
#############################################################################################################
hotspot.scan <- function(cross, scan.hl, lod.thr, quant.level, window = 5, verbose = FALSE)
{
  ## Kludge to get chr and pos[
  cross$pheno$trait <- rnorm(nind(cross))
  if(is.null(cross$geno[[1]]$prob))
    cross <- calc.genoprob(cross, step=0.5, err = 0.002,
                           map.function = "c-f", stepwidth = "max")
  
  chr.pos <- scanone(cross, pheno.col = find.pheno(cross, "trait"))[,1:2]

  hot.scan <- pull.hotspots(cross, scan.hl, chr.pos, lod.thr, quant.level)

  if(verbose) cat("make.maxlod\n")
  max.hl <- make.maxlod(scan.hl, chr.pos)

  if(verbose) cat(paste("smoth.neqtl", window, "\n"))
  neqtl.pos <- smooth.neqtl(scan.hl, chr.pos, max.hl, lod.thr, window)
  hot.scan$max.N.window <- neqtl.pos$nqtl

  if(verbose) cat(paste("smoth.neqtl", 1, "\n"))
  neqtl.pos0 <- smooth.neqtl(scan.hl, chr.pos, max.hl, lod.thr, 1)
  hot.scan$max.N0 <- neqtl.pos0$nqtl

  hot.scan
}
scan.hl.plot <- function(scan.hl)
{
  for(i in levels(scan.hl$chr)) {
    print(xyplot(lod~pos, scan.hl[scan.hl$chr == i,], group = phenos, type = "l",
                 main = paste("chromosome", i)))
  }
  invisible()
}
hotspot.plot <- function(hot.scan, quant.thr = NULL, maps = NULL, main = "")
{
  for(i in levels(hot.scan$chr)) {
    ylims <- range(c(hot.scan[hot.scan$chr == i, -(1:2)]))
    plot(hot.scan, lodcolumn=1, chr = i, col = "black", ylab = "hotspot size", ylim = ylims)
    mtext(paste(main, "raw (black) smoothed (red) sliding (blue)"), 3, 3)
    plot(hot.scan, lodcolumn=3, chr = i, col = "red", add = TRUE)
    plot(hot.scan, lodcolumn=2, chr = i, col = "blue", add = TRUE)

    if(!is.null(maps))
      add.rug(i, "", maps, use.cM = TRUE)

    if(!is.null(quant.thr))
      abline(h = quant.thr, lwd = 2, lty = 2, col = "red")
  }
  
  invisible(hot.scan)
}
