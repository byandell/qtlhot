## see ~/p/private/diabetes1/diabetes10/scan.perm/func1.R
lininterp <- function(x, y, xnew, ynew)
{
  if(!missing(ynew)) {
    if(!missing(xnew)) 
      stop("Give one of xnew or ynew, not both.")
    return(lininterp(y, x, ynew))
  }

  if(missing(xnew))
    stop("Give one of xnew or ynew.")

  wh <- is.nan(x) | is.nan(y)
  if(any(wh)) {
    x <- x[!wh]
    y <- y[!wh]
  }
  if(length(x) != length(y))
    stop("x and y must have the same length")
  if(length(x) < 2)
    stop("x and y must have length > 1")
  if(any(diff(x) < 0) || any(diff(y) < 0))
    stop("x and y must be non-decreasing")
  
  xd <- diff(range(x))
  yd <- diff(range(y))

  z <- rep(NA, length(xnew))
  seen <- rep(FALSE, length(xnew))

  # on top
  wh <- match(xnew, x)
  if(any(!is.nan(wh))) {
    z[!is.nan(wh)] <- y[wh[!is.nan(wh)]]
    seen[!is.nan(wh)] <- TRUE
  }
    
  # to the right
  wh <- !is.nan(xnew) & xnew > max(x) & !seen
  if(any(wh)) {
    z[wh] <- (xnew[wh]-max(x))/xd*yd+max(y)
    seen[wh] <- TRUE
  }
  
  # to the left
  wh <- !is.nan(xnew) & xnew < min(x) & !seen
  if(any(wh)) {
    z[wh] <- min(y) - (min(x)-xnew[wh])/xd*yd
    seen[wh] <- TRUE
  }

  # in between
  wh <- !is.nan(xnew) & xnew >= min(x) & xnew <= max(x) & !seen
  if(any(wh)) {
    seen[wh] <- TRUE

    for(i in which(wh)) {
      xl <- max(x[x <= xnew[i]])
      xr <- min(x[x >= xnew[i]])
      yl <- max(y[x <= xnew[i]])
      yr <- min(y[x >= xnew[i]])
      z[i] <- (xnew[i]-xl)/(xr-xl)*(yr-yl) + yl
    }
  }
  
  z
}
  
smoothall <- function(themax, thechr, thepos, window=5)
{  
  thesmooth <- vector("list", length(themax))
  names(thesmooth) <- names(themax)
  for(i in names(themax))
      thesmooth[[i]] <- smoothchr(themax[[i]], thepos[thechr==i], window=window)
   out <- NULL
  for(i in 1:length(thesmooth))
    out <- rbind(out, data.frame(chr=rep(names(themax)[i], nrow(thesmooth[[i]])),
                    pos=thesmooth[[i]][,1], nqtl=thesmooth[[i]][,2]))
  class(out) <- c("scanone", "data.frame")

  ## This way loses the marker information.
  rownames(out) <- paste("c", out[,1], ".loc", 1:nrow(out), sep="")

  out
}

## Uses positions from thepos for smoothing: ATB 9/10/09 ##
smoothchr <- function(themax, thepos, window=5)
{
  #theloc <- sort(unique(c(thepos, seq(0, max(thepos)), by=0.2)))
  theloc <- sort(thepos)

  temploc <- c(themax, theloc)
  tempval <- c(rep(1, length(themax)), rep(0, length(theloc)))
  o <- order(temploc)
  temploc <- temploc[o]
  tempval <- tempval[o]
#  smoothed <- runningmean(temploc, tempval, at=theloc, window=window, what="sum") 
  smoothed <- runningmean(temploc, tempval, window=window, what="sum")
  u <- tempval == 0
  return(cbind(pos=temploc[u], smoothed[u]))
}

findpeaks <- function(results, lodcolumn=7, window=5, n.peaks=20)
{
  results$chr <- as.character(results$chr)
  for(i in 1:n.peaks) {
    temp <- max(results, lodcolumn=lodcolumn)
    if(i==1) 
      output <- temp
    else
      output <- rbind(output, temp)
    
    chr <- as.character(temp[[1]])
    pos <- temp[[2]]
    results[results$chr==chr & results$pos >= pos-window & results$pos <= pos+window,lodcolumn+2] <- NA
  }
  output
}
      
