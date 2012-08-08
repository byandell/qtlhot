##############################################################################
OGetCisCandReg <- function(cand.reg, scan, lod.thr, drop = 1.5) {
  all.trait.nms <- cand.reg[, 1]
  cand.reg <- cand.reg[cand.reg[, 2] == cand.reg[, 4],]
  n <- nrow(cand.reg)
  trait.nms <- names(scan)[-c(1, 2)]
  is.cis <- rep(FALSE, n)
  peak.pos.lower <- rep(NA, n)
  peak.pos.upper <- rep(NA, n)
  for (i in 1 : n) {
    phys.pos <- cand.reg[i, 3]
    trait.index <- which(trait.nms == cand.reg[i, 1]) + 2
    sscan <- scan[scan[, 1] == cand.reg[i, 4], c(1, 2, trait.index)]
    lod.interval <- lodint(sscan, chr = cand.reg[i, 4], drop = drop)
    peak.pos.lower[i] <- lod.interval[1, 2]
    peak.pos.upper[i] <- lod.interval[nrow(lod.interval), 2]
    lod.phys.pos <- sscan[which.min(abs(sscan[, 2] - cand.reg[i, 3])), 3]
    if (phys.pos >= peak.pos.lower[i] & phys.pos <= peak.pos.upper[i] & 
        lod.phys.pos >= lod.thr) {
      is.cis[i] <- TRUE
    }
  }
  out <- data.frame(cand.reg, peak.pos.lower, peak.pos.upper)
  index <- match(out[is.cis, 1], all.trait.nms)
  list(cis.reg = out[is.cis,], cis.index = index)
}
##############################################################################
OGetCandReg <- function(scan, annot, traits, lod.thr, drop = 1.5,
                        nms = names(scan)[-c(1,2)])
  
{
  ## want to use highlod instead of scan below.
  ## need to decode highlod.
  ## currently this only gets max over genome; want max by chr, yes?
  traits <- unique(traits)
  n <- length(traits)
  out <- data.frame(matrix(NA, n, 6))
  names(out) <- c("gene", "phys.chr", "phys.pos", "peak.chr", "peak.pos",
                  "peak.lod")

  for (i in 1:n) {
    ii <- match(traits[i], annot[, 1])
    out[i, 1:3] <- annot[ii, c(1, 3, 5)]
    trait.chr <- annot[ii, 3]
    trait.pos <- annot[ii, 5]
    if (!is.na(trait.pos)) {
      peak.index <- which.max(scan[, traits[i]])
      peak.lod <- scan[peak.index, traits[i]]
      peak.chr <- scan[peak.index, 1]
      peak.pos <- scan[peak.index, 2]
      if (peak.lod >= lod.thr) {
        out[i, 4] <- peak.chr
        out[i, 5] <- peak.pos
        out[i, 6] <- peak.lod
      }
    }     
    cat(" ", i, "\n")
  }
  subset(out, !is.na(out[, 4]))
}
##############################################################################
OGetCoMappingTraits <- function(traits, scan, lod.thr, drop = 1.5) {
  n <- nrow(traits)
  out <- vector(mode = "list", length = n)
  names(out) <- traits[, 1]
  nms1 <- names(scan)[-c(1, 2)]
  for (i in 1:n) {
    peak.pos <- traits[i, 5]
    sub.scan <- scan[scan[, 1] == traits[i, 4], ]
    peak.lods <- apply(sub.scan[, -c(1, 2)], 2, max)
    aux <- which(peak.lods >= lod.thr)
    aux <- aux[-match(traits[i, 1], names(aux))]
    nms2 <- names(aux)
    is.target <- rep(FALSE, length(nms2))
    for (j in 1:length(nms2)) {
      trait.index <- match(nms2[j], nms1)
      sscan <- sub.scan[, c(1, 2, trait.index + 2)]
      lod.interval <- lodint(sscan, chr = traits[i, 4], drop = drop)
      lb <- lod.interval[1, 2]
      ub <- lod.interval[nrow(lod.interval), 2] ## was 3; error if there are tied LOD scores.
      if(peak.pos >= lb & peak.pos <= ub)
        is.target[j] <- TRUE
    }
    out[[i]] <- nms2[is.target]
    cat(" ", i, "\n")
  }
  out
}
