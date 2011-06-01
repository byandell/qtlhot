## Generates a "null dataset" cross
##
sim.null.cross <- function(chr.len = rep(400,16), n.mar=185, n.ind = 112, type = "bc",
                           n.pheno = 6000, latent.eff = 1.5, res.var = 1,
                           init.seed = 92387475)
{
  set.seed(init.seed)
  mymap <- sim.map(len = chr.len, n.mar = n.mar, include.x=FALSE, eq.spacing=TRUE)
  scross <- sim.cross(map = mymap, n.ind = n.ind, type = type)
  scross <- calc.genoprob(scross, step=0)

  sim.null.pheno.data(cross = scross, n.pheno = n.pheno, latent.eff = latent.eff, res.var = res.var)
}
sim.null.pheno.data <- function(cross, n.pheno, latent.eff, res.var)
{
  n <- nind(cross)
  latent <- rnorm(n, 0, sqrt(res.var))
  ErrorM <- matrix(rnorm(n * n.pheno, 0, sqrt(res.var)), n, n.pheno)
  pheno <- data.frame(latent*latent.eff + ErrorM)
  names(pheno) <- paste("P", 1:n.pheno, sep="")
  cross$pheno <- pheno
  
  cross
}

## Generate hotspots for the simulated examples in the manuscript.
include.hotspots <- function(cross,
                             hchr,
                             hpos,
                             hsize,
                             Q.eff,
                             latent.eff,
                             lod.range.1,
                             lod.range.2,
                             lod.range.3,
                             res.var=1,
                             nT,
                             init.seed)
{
  get.closest.pos.nms <- function(pos, cross, chr) {
    ## map <- attributes(cross$geno[[chr]]$prob)$map
    ## This can be simplified.
    map <- pull.map(cross, chr)
    q.nms <- names(map)
    map.pos <- as.numeric(map)
    tmp <- which.min(abs(map.pos-pos))
    closest.pos <- map.pos[tmp]
    nms <- q.nms[tmp]
    list(nms,closest.pos)
  }
  pull.prob <- function(cross) {
    out <- vector(mode = "list", length = nchr(cross))
    names(out) <- names(cross$geno)
    for(i in names(out))
      out[[i]] <- cross$geno[[i]]$prob  
    out
  }
  
  set.seed(init.seed)
  
  hk.prob <- pull.prob(cross)

  update.pheno <- function(cross, hchr, hpos, hsize, Q.eff, latent.eff, lod.range,
                           res.var, index, hk.prob) {
    map <- attributes(hk.prob[[hchr]])$map
    q.nms <- names(map)
    map.pos <- as.numeric(map)
    tmp <- which.min(abs(map.pos-pos))
    closest.pos <- map.pos[tmp]
    N.pos <- q.nms[tmp]
    M.pos <- get.closest.pos.nms(hpos, cross, hchr)[[1]]
    
    M.dummy <- hk.prob[[hchr]][, M.pos, 1] - hk.prob[[hchr]][, M.pos, 2]
    M <- M.dummy * Q.eff + rnorm(nind(cross), 0, sqrt(res.var))

    ## QTL effect
    lod <- runif(hsize, lod.range[1], lod.range[2])
    r2 <- 1 - 10 ^ (-2 * lod / nind(cross))
    beta <- sqrt(r2 * res.var * (1 + latent.eff ^ 2) / (Q.eff ^ 2 * (1 - r2) - res.var * r2))

    for(j in seq(length(index)))
      cross$pheno[, index[j]] <- beta[j] * M + cross$pheno[, index[j]]

    cross
  }

  ## Why 50 for first, 500 for 2nd and 3rd?
  ## Why strange lod.range for 2nd?
  index1 <- sample(1:nT, 50, replace = FALSE)
  cross <- update.pheno(cross, hchr[1], hpos[1], hsize[1], Q.eff, latent.eff,
                        lod.range.1, res.var, index1, hk.prob)

  index2 <- sample((1:nT)[-index1], 500, replace = FALSE)
  cross <- update.pheno(cross, hchr[2], hpos[2], hsize[2], Q.eff, latent.eff,
                        c(lod.range.1[2], lod.range.2[1]), res.var, index2, hk.prob)

  index3 <- sample((1:nT)[-c(index1, index2)], 500, replace = FALSE)
  cross <- update.pheno(cross, hchr[3], hpos[3], hsize[3], Q.eff, latent.eff,
                        lod.range.3, res.var, index3, hk.prob)

  cross
}
