## Directory path is indicated with "dirpath" character string
##
## I/O files are in internal R format (RData) using load/save commands
## except "groups.txt", which is used to set width of parallelization.
##
## I/O files are on AFS in
## /u/y/a/yandell/public/html/sysgen/qtlhot/condor
## or equivalently via URL
## http://www.stat.wisc.edu/~yandell/sysgen/qtlhot/condor
##
####################################################################################
parallel.qtlhot <- function(phase, index = 1, ..., dirpath = ".")
{
  switch(phase,
         qtlhot.phase1(dirpath, index, ...),
         qtlhot.phase2(dirpath, index, ...),
         qtlhot.phase3(dirpath, index, ...),
         parallel.error(1, phase, index))
  
  parallel.error(0, phase, index)
}
parallel.error <- function(num, phase = 0, index = 1)
{
  ## See file errorcodes.txt for explanation.
  if(phase == 5)
    outname <- "RESULT"
  else
    outname <- paste("RESULT", phase, index, sep = ".")
  write.table(num, file = outname,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  if(num)
    stop(parallel.message(num))

  invisible()
}
parallel.message <- function(num)
{
  msg <- as.matrix(c("OK",
                     "phase number not between 1 and 3",
                     "index to groups must be supplied as line number for file groups.txt",
                     "index to groups must be integer (line number for file groups.txt)",
                     "index to groups must be between 1 and number of lines in groups.txt",
                     "no bic RData files found",
                     "index to MCMC runs must be supplied",
                     "index to groups must be integer",
                     "index to groups must be between 1 and nruns",
                     "no MCMC RData files found",
                     "number of MCMC runs does not match nruns",
                     "number of BIC runs does not match groups"))
  dimnames(msg)[[1]] <- seq(0, nrow(msg) - 1)
  if(missing(num))
    write.table(msg, file = "errorcodes.txt",
                quote = FALSE, col.names = FALSE)
  else
    msg[as.character(num), 1]
}
####################################################################################
qtlhot.phase0 <- function(dirpath, init.seed = 92387475,
                          len = rep(400,16),
                          n.mar = 185, n.ind = 112,
                          n.phe = 100, latent.eff = 0, res.var = 1,
                          lod.thrs = c(4.63,4.17,3.93,3.76,3.65,3.56,3.47,3.39,3.34,3.30),
                          ...)
{
  ## PHASE 0: Cross object initialization. Create file "cross.RData" for Phase1,Phase3.
  
  set.seed(init.seed)

  mymap <- sim.map(len = len, n.mar = n.mar, include.x = FALSE, eq.spacing = TRUE)
  cross <- sim.cross(map = mymap, n.ind = n.ind, type = "bc")
  cross <- calc.genoprob(cross, step=0)
  
  ## Simulate new phenotypes.
  cross <- sim.null.pheno.data(cross, n.phe, latent.eff, res.var)
  
  save(cross, lod.thrs, file = file.path(dirpath, "cross.RData"), compress = TRUE)
}
####################################################################################
qtlhot.phase1 <- function(dirpath, index = 0,
                          params.file = "params.txt",
                          cross.file = "cross.RData",
                          cross.name = "cross",
                          n.ind = nind(cross),
                          lod.thrs = c(4.63,4.17,3.93,3.76,3.65,3.56,3.47,3.39,3.34,3.30),
                          alpha.levels = c(1:10) / 100, ## Nominal significance levels.
                          Nmax = 100, ## Maximum hotpsot size recorded.
                          n.perm = 1000, ## Number of permutations.
                          n.split = 100, ## Number of splits of permutations.
                          latent.eff = 1.5, ## Latent effect determines correlation among traits.
                          res.var = 1, ## Residual variance.
                          n.phe = nphe(cross), ## Number of traits.
                          nruns = 1,
                          big = FALSE,
                          ...)
{
  ## PHASE 1: Set up cross object. Needed in phases 2.
  ##
  ## Input files:
  ##       cross.RData: cross object and possibly lod.thrs or other objects.
  ##       params.txt
  ##
  ## Output files:
  ##       Phase1.RData
  ##

  ## Get any parameters in file. These will overwrite passed arguments.
  eval(parse(file.path(dirpath, params.file)))

  ## Cross object. Load if not done already.
  if(!exists(cross.name))
    load(file.path(dirpath, cross.file))
  
  ## Change name of cross object to "cross" for internal use.
  if(cross.name != "cross")
    cross <- get(cross.name)

  ## cross.index is used when multiple phase1 jobs are spawned.
  cross.index <- as.numeric(index)

  if(big) {
    ## Used for big data run.
    ## Each run is separate permutation.
    ## Make sure n.split is equal to n.perm. n.perm set to 1 in big.phase1.
    n.split <- max(1, n.perm)
    big.phase1(dirpath, cross.index, params.file, cross, lod.thrs, Nmax, n.perm,
               n.split, ...)
  }
  else {
    ## Used for studying properties of qtlhot.

    n.perm <- ceiling(n.perm / n.split)
      
    if(cross.index > 0 & nruns > 1) {
      ## Keep markers but re-simulate genotypes each time.
      mymap <- pull.map(cross)
      cross <- sim.cross(map = mymap, n.ind = n.ind, type = class(cross)[1])
      cross <- calc.genoprob(cross, step=0)
      
      ## Simulate new phenotypes.
      cross <- sim.null.pheno.data(cross, n.phe, latent.eff, res.var)
    }
    else
      cross.index <- 0

    Nmax <- min(Nmax, n.phe)
    
    ## Save all relevant objects for later phases.
    save(cross, n.phe, latent.eff, res.var, Nmax, n.perm, n.split,
         alpha.levels, lod.thrs, cross.index, big,
         file = file.path(dirpath, "Phase1.RData"),
         compress = TRUE)
  }

  ## Need to write a file with n.groups lines and group.size columns.
  ## The file groups.txt is examined by parallelizer (SOAR).
  write.table(c(cross.index, n.split),
              file = file.path(dirpath, "groups.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}
####################################################################################
qtlhot.phase2 <- function(dirpath, index = NULL, ..., big = FALSE, verbose = FALSE)
{
  ## PHASE 2: NL,N and WW permutations. 
  ##          Slow. Run on condor nodes. Sized by n.perm.
  ##
  ## Steps for time calculations: n.perm * nind(cross) * nmar(cross) * nphe(cross)
  ## Input files:
  ##       Phase1.RData: cross, n.phe, n.perm, lod.thrs, 
  ##
  ## Output file (one per invocation):
  ##       permi.RData
  ##

  cross.index <- scan("groups.txt", 0)[1]

  ## Load Phase 1 computations.
  infile <- "Phase1.RData"
  load(file.path(dirpath, infile))

  ## Quality check of index.
  if(missing(index))
    parallel.error(2, 2, index)
  index <- as.integer(index)
  if(is.na(index))
    parallel.error(3, 2, index)
  if(index < 1 | index > n.split)
    parallel.error(4, 2, index)

  if(big)
    return(big.phase2(dirpath, index))

  outfile <- paste("perm", ".", cross.index, "_", index, ".RData", sep = "")

  ## Following is in NL.N.perm stuff.
  ## sum.threshold is big loop. Have perm loop within it. n.phe*n.perm runs of scanone per dataset.
  ## For simulation study, have many datasets!

  ## Creates max.N of size n.perm x n.lod and max.lod.quant of size n.perm x Nmax.
  ## Size of n.perm determines the run time.
  NLN <- NL.N.perm(cross, Nmax, n.perm, lod.thrs, drop = drop,
                   verbose = verbose)

  mycat("scanone", verbose)
  scanmat <- scanone(cross, pheno.col = seq(n.phe), method = "hk")
  chr <- scanmat[,1]
  scanmat <- as.matrix(scanmat[, -(1:2), drop = FALSE])
  
  ## we use the lowest lod threshold in the LOD drop interval 
  mycat("set.to.zero.beyond.drop.int", verbose)
  scan.drop <- set.to.zero.beyond.drop.int(chr, scanmat, min(lod.thrs), drop)
  
  WW <- WW.perm(scan.drop, lod.thrs, n.perm)
  
  save(NLN, WW, index, lod.thrs, alpha.levels, drop,
       file = outfile, compress = TRUE)
}
####################################################################################
qtlhot.phase3 <- function(dirpath, index = NULL, ...,
                          dirpath.save = dirpath, big = FALSE, verbose = FALSE)
{
  ## PHASE 3: Sample Markov chain (MCMC). Parallelize.
  ##          Fast: Run on scheduler.
  ##          Slow. Run on condor nodes.
  ##
  ## Input files:
  ##       Phase1.RData
  ##       permN.RData (multiple files with N = 1,...,n.split
  ##       cross.RData: cross object and possibly lod.thrs or other objects.
  ##
  ## Output files (one per invocation):
  ##       Phase3.RData: max.N, max.lod.quant
  ##
  ## See Phase 2 for explanation of NLNpermi.RData files.
  ## All NLNpermi.RData files are combined to make max.N, max.lod.quant.

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  if(big)
    return(big.phase3(dirpath, index, cross.index))

  ## This could be done once, but it would require splitting this phase in two.
  ## Besides, it is quite fast.
  ## Read in saved BIC scores and combine into one object.
  outfile <- paste("perm", ".", cross.index, "_", ".*RData", sep = "")
  filenames <- list.files(dirpath, outfile)
  if(!length(filenames))
    parallel.error(5, 3, index)
  if(length(filenames) != n.split)
    parallel.error(11, 3, index)

  n.perms <- n.perm * n.split
  n.lod <- length(lod.thrs)
  
  max.N <- max.WW <- matrix(NA, n.perms, n.lod)
  max.lod.quant <- matrix(NA, n.perms, Nmax)

  i.perm <- seq(n.perm)
  for(i in seq(length(filenames))) {
    load(file.path(dirpath, filenames[i]))

    ## Do any quality checking here.
    
    max.N[i.perm, ] <- NLN$max.N
    max.lod.quant[i.perm, ] <- NLN$max.lod.quant
    max.WW[i.perm, ] <- WW
    
    i.perm <- i.perm + n.perm
  }

  NL.N.thrs <- NL.N.summary(max.lod.quant, max.N, alpha.levels)
  WW.thrs <- WW.summary(max.WW, alpha.levels)

  ## Now compare permutations to original cross object from Phase1.

  out.sim <- sum.threshold(cross, n.phe, latent.eff, res.var,
                           lod.thrs, drop, seq(Nmax), n.perms, alpha.levels,
                           NL.N.thrs, WW.thrs,
                           verbose)

  outfile <- paste("Phase3", ifelse(cross.index > 0, cross.index, ""), ".RData", sep = "")
  
  save(latent.eff, Nmax, out.sim, 
       file = file.path(dirpath.save, outfile),
       compress = TRUE)
}
