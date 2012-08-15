######################################################################
# hotsize.R
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
# Contains: hotsize, hotsize.scanone, hotsize.highlod,
#           print.hotsize, summary.hotsize, plot.hotsize
######################################################################
hotsize <- function(hotobject, ...) UseMethod("hotsize")

hotsize.scanone <- function(hotobject, lod.thr = NULL, drop.lod = 1.5, ...)
{
  hotsize.highlod(highlod(hotobject, lod.thr, drop.lod), lod.thr, ...)
}
hotsize.highlod <- function(hotobject, lod.thr = NULL, window = NULL, quant.level = NULL, ...)
{
  if(length(lod.thr) > 1)
    stop("hotsize only allows one lod.thr value")
  
  scan <- hotobject$chr.pos
  
  highlod <- hotobject$highlod
  if(is.null(lod.thr))
    lod.thr <- attr(hotobject, "lod.thr")
  
  if(!is.null(lod.thr)) {
    highlod <- highlod[highlod$lod >= lod.thr,, drop = FALSE]
    attr(scan, "lod.thr") <- lod.thr
  }
  
  if(!nrow(highlod))
    return(NULL)
  
  ## Straight count of LODs above threshold. Includes shoulders and peaks.
  tbl <- table(highlod$row)
  tmp <- rep(0, nrow(scan))
  if(length(tbl))
    tmp[as.numeric(names(tbl))] <- tbl
  scan$max.N <- tmp
  
  ## Smoothed count of peak LODs per chr only above threshold, smoothed by window.
  if(!is.null(window)) {
    window <- round(window)
    if(window > 0) {
      scan$max.N.window <- smooth.neqtl(highlod, chr.pos, window = window)$nqtl
      attr(scan, "window") <- window
    }
  }
  
  if(!is.null(quant.level)) {
    ## If sliding thresholds supplied.
    if(!is.null(lod.thr))
      quant.level <- quant.level[quant.level >= lod.thr & !is.na(quant.level)]
    
    if(length(quant.level)) {
      ## Want to work down quant.level. One way is to see how many are above smallest
      ## and stop if above the index.
      index <- rep(TRUE, nrow(highlod))
      hot.size <- rep(0, nrow(scan))
      for(hot.crit in rev(seq(length(quant.level)))) {
        if(verbose)
          cat(hot.crit, "\n")
        above <- highlod$lod[index] >= quant.level[hot.crit]
        if(!any(above))
          break;
        tbl <- table(highlod$row[index][above])
        tbl <- tbl[tbl >= hot.crit]
        if(length(tbl)) {
          ## Record hotspot size if above hot.crit for thr, then mask those loci out.
          hot.size[as.numeric(names(tbl))] <- tbl
          index[highlod$row %in% names(tbl)] <- FALSE
        }
        if(!any(index))
          break;
      }
      scan$quant <- hot.size
      attr(scan, "quant.level") <- quant.level
    }
  }
  class(scan) <- c("hotsize", class(scan))
  scan
}
#############################################################################################
print.hotsize <- function(object, ...) print(summary(object, ...))
summary.hotsize <- function(object, ...)
{
  
  cat("hotsize elements: ", paste(names(object)), "\n")

  lod.thr <- attr(object, "lod.thr")
  if(!is.null(lod.thr))
    cat("LOD threshold:", lod.thr, "\n")
  window <- attr(object, "window")
  if(!is.null(window))
    cat("smooth window:", window, "\n")
  quant.level <- attr(object, "quant.level")
  if(!is.null(quant.level)) {
    cat("quantile level summary:\n")
    summary(quant.level)
  }
  cat("\n")
  format <- ifelse(ncol(object ==3), "onepheno", "allpeaks")
  NextMethod(object, format = format, ...)
}    
#############################################################################################
plot.hotsize <- function(x, ylab = "counts", quant.axis = pretty(c(1, length(quant.level))),
             col = 1:3, ...)
{
  ## Use NextMethod, but repeatedly.
  class(x) <- class(x)[-1]
  ## max.N
  plot(x, lodcolumn = 1, ylab = ylab, col = col[1], ...)

  ## max.N.window
  window <- attr(x, "window")
  if(!is.null(window)) {
    plot(x, lodcolumn = 2, col = col[2], ...)
  }

  ## quant
  quant.level <- attr(x, "quant.level")
  if(!is.null(quant.level)) {
    plot(x, lodcolumn = match("quant", names(x)) - 2, col = col[3],
                 add = TRUE, ...)

    ## Add right axis for quantile LOD level.
    if(length(quant.axis)) {
      axis(4, at = quant.axis, label = round(quant.level[quant.axis], 1), las = 1, cex = 0.35)
      mtext("sliding LOD thresholds", 4, 1, cex = 1.5)
    }
  }
  invisible()
}
