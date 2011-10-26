## Routines to show hotspots. None here yet.
pull.hotspots <- function(cross, scan.hl, chr.pos = NULL, lod.thr = 5, slide.thr = NULL)
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
    chr.pos
    if(!is.null(slide.thr)) {
      ## If sliding thresholds supplied.
      slide.thr <- slide.thr[slide.thr >= lod.thr & !is.na(slide.thr)]
      if(length(slide.thr)) {
        ## Want to work down slide.thr. One way is to see how many are above smallest
        ## and stop if above the index.
        index <- rep(TRUE, nrow(scan.hl))
        hot.size <- rep(0, nrow(chr.pos))
        for(hot.crit in rev(seq(length(slide.thr)))) {
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
