## This function takes a scanone object (scan) and for each trait and each
## chromosome it: (1) finds out the max LOD peak; (2) if the max LOD value
## is smaller than the map threshold (thr) it set to zero all LOD values;
## (3) if the max LOD value is higher than thr, it computes the LOD drop 
## interval around the peak and set to zero all LOD values outside the LOD 
## LOD interval.
## 
set.to.zero.beyond.drop.int <- function(chr, scanmat, thr, drop = 1.5)
{
  mylodint <- function(x, drop = 1.5)
  {
    drops <- x < (max(x) - drop)
    if(any(drops))
      x[drops] <- 0
    x
  }

  mychr <- levels(chr)
  for(i in mychr){
    pointer <- which(scanmat[,1] == i)
    if(length(pointer)) {
      if(max(scanmat[pointer,]) < thr)
        ## Zero out if no peak above threshold.
        scanmat[pointer,] <- 0
      else
        ## Zero out below drop from peak.
        scanmat[pointer,] <-
          apply(scanmat[pointer,, drop = FALSE], 2, mylodint, drop)
    }
  }
  scanmat
}
