######################################################################
# diag.R
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
# Contains: myplot.err, sliding.bar.plot
######################################################################

myplot.err <- function(out.err, dim.num = 1, ylab = deparse(substitute(out.err)), jitter = FALSE)
{
  mylevels <- lapply(dimnames(out.err), as.numeric)
  mylen <- sapply(mylevels, length)
  dim.two <- 3 - dim.num
  
  if(dim.num == 1) {
    myfun <- function(i) i + mylen[1] * seq(0, mylen[2] - 1)
    xlab <- "lod threshold by alpha level"
  }
  else {
    myfun <- function(i) seq(1, mylen[1]) + mylen[1] * (i - 1)
    xlab <- "alpha level by lod threshold"
  }

  plot(range(mylevels[[dim.num]]), range(c(out.err)), type = "n",
       xlab = xlab, ylab = ylab)
  for(i in seq(mylen[dim.num]))
    lines(mylevels[[dim.num]],
          if(jitter) jitter(out.err[myfun(i)]) else out.err[myfun(i)],
          type = "b", pch = as.character(i))
}

## Generates the sliding bar figure (not an image anymore)
##
sliding.bar.plot <- function(scan, lod.thr, size.thr, gap=50, y.axes=NULL)
{
  ###
  sliding.bar.matrix <- function(scan, lod.thr, size.thr)
  {
    counts <- count.thr(scan, lod.thr, TRUE)
    out <- matrix(FALSE,nrow(counts),ncol(counts))
    for(i in order(size.thr))
      out[ i, counts[i,] > size.thr[i] ] <- TRUE
    out
  }
  ###
  matrix.to.plot <- function(M, map, myrug)
  {
    pM <- matrix(FALSE, nrow(M), ceiling(max(myrug)))
    rrug <- as.numeric(round(myrug,0))
    for(i in 1:ncol(M))
      pM[,rrug[i]] <- M[,i]
    pM
  }
  ###
  create.rug <- function(map, gap)
  {
    chrs <- unique(map[,1])
    nchrs <- length(chrs)
    maxpos <- rep(NA,nchrs)
    myrug <- map[map[,1]==1,2]
    chr.legend.pos <- median(myrug)
    for(i in 2:length(chrs)){
      aux <- max(myrug) + gap
      myrug <- c(myrug, map[map[,1]==i,2] + aux)
      chr.legend.pos <- c(chr.legend.pos, median(map[map[,1]==i,2] + aux))
    }
    list(myrug=myrug, chr.legend.pos=chr.legend.pos)
  }
  ###
  N <- max(size.thr)
  sbm <- sliding.bar.matrix(scan, lod.thr, size.thr)
  map <- scan[,1:2]
  myrug <- create.rug(map, gap)
  M <- matrix.to.plot(sbm, map, myrug[[1]])
  lod.thr <- as.character(round(lod.thr,2))
  xaxis <- c(1:ncol(M))
  par(mar=c(5, 4, 4, 5) + 0.1) 
  plot(xaxis, xaxis, type="n", ylim=c(0,N), xaxt="n", xlab="Chromosome",
       yaxt="n", ylab="Hotspot size", cex.lab=1.5)
  rug(myrug[[1]], 0.02, quiet=TRUE)
  axis(side=1, labels=as.character(unique(map[,1])), at=myrug[[2]], 
       cex.axis=1.5, tick=FALSE)
  if(is.null(y.axes))
    y.axes <- quantile(c(1:N), c(0, 0.25, 0.50, 0.75, 1))
  axis(side=2, labels=as.character(y.axes), at=y.axes, cex.axis=0.9, las=1)
  axis(side=4, labels=lod.thr[y.axes], at=y.axes, cex.axis=0.9, las=1)
  mtext("LOD threshold",side=4,cex=1.4,line=3.5,adj=0.545)
  for(i in 1:N){
    for(j in 1:ncol(M)){
      if(M[i,j]) segments(x0=j, x1=j, y0=i-1/2, y1=i+1/2, lwd=0.1)
    }
  }
}

