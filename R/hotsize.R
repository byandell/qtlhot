hotsize <- function(object, ...) UseMethod("hotsize")

pull.hotspots <- function(cross, scan.hl, chr.pos = NULL, lod.thr = 5, slide.thr = NULL, verbose = FALSE)
hotsize.highlod <- function(object, lod.thr = NULL, quant.level = NULL, ...)
{
  scan <- object$chr.pos

  highlod <- object$highlod
  if(!is.null(lod.thr))
    highlod <- highlod[highlod$lod >= lod.thr,, drop = FALSE]

  if(nrow(highlod)) {
    ## Straight count of LODs above threshold. Includes shoulders and peaks.
    tbl <- table(highlod$row)
    tmp <- rep(0, nrow(scan))
    if(length(tbl))
      tmp[as.numeric(names(tbl))] <- tbl
    scan$max.N <- tmp

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
  else
    NULL
}
#############################################################################################
plot.hotsize(x, ylab = "counts", quant.axis = pretty(c(1, length(quant.level))), ...)
{
  ## NextMethod:
  plot.scanone(x, lodcolumn = 1, ylab = ylab, ...)

  ## Add right axis for quantile LOD level.
  quant.level <- attr(hot.scan1, "quant.level")
  if(!is.null(quant.level)) {
    ## Need to check lod columns here!
    plot.scanone(x, ylab = ylab, lodcolumn = 2, add = TRUE, ...)
    if(length(quant.axis)) {
      axis(4, at = quant.axis, label = round(quant.level[quant.axis], 1), las = 1, cex = 0.35)
      mtext("sliding LOD thresholds", 4, 1, cex = 1.5)
    }
  }
  invisible()
}
