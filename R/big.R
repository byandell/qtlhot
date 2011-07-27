big.phase1 <- function(dirpath = ".", cross.index = 0, params.file,
                       nruns, cross, lod.thrs, Nmax = 2000, n.perm = 1,
                       trait.file, trait.matrix, trait.index, droptrait.names, batch.effect = NULL,
                       drop.lod = 1.5, lod.min = min(lod.thrs),
                       size.set = 250, window = 5, seed = 0, ...,
                       verbose = FALSE)
{
  ## Want this to depend on parallel.phase1 as much as possible.
  ## Just include new stuff here.
  
  ## Phase 1 for big dataset.
  ## Sets up initial permutations.
  
  ## Get any parameters in file. These will overwrite passed arguments.
  eval(parse(file.path(dirpath, params.file)))

  ## Calculate genotype probabilities.
  cross <- calc.genoprob(cross, step=0.5, err = 0.002,
                         map.function = "c-f", stepwidth = "max")
  
  ## Subset to individuals with batch if used.
  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))
  
  ## Get trait values for selected traits.

  attach(file.path(dirpath, trait.file))
  trait.names <- dimnames(get(trait.matrix))[[2]]
  detach()
  all.traits <- seq(trait.names)[-match(droptrait.names, trait.names)]

  cross.index <- as.numeric(cross.index)

  ## Random numbers.
  if(n.perm > 1) {
    if(seed[1] > 0)
      set.seed(seed[1])

    n.ind <- nind(cross)
    seeds <- sample(c(98765:987654321), n.perm, replace = FALSE)
  }
  else
    seeds <- rep(0, n.perm)

  save(cross, n.perm, seeds, Nmax, drop.lod, lod.thrs, lod.min, batch.effect,
       control.probes, cross, trait.file, trait.matrix, trait.index,
       trait.names, all.traits, size.set, window, cross.index,
       file = "Phase1.RData", compress = TRUE)
  invisible()
}
#############################################################################################
big.phase2 <- function(dirpath = ".", index, ..., verbose = FALSE)
{
  ## Phase 2.
  ## Loop on sets of phenotypes, adding only what is needed.
  ## Loop on permutations internal to scanone.permutations.
  
  load("Phase1.RData")
  
  covars <- sexbatch.covar(cross, batch.effect)

  n.ind <- nind(cross)
  n.phe <- nphe(cross)
  n.traits <- length(all.traits)
  num.sets <- ceiling(n.traits / size.set)
  trait.nums <- seq(size.set)

  if(n.perm > 1) {
    ## Random permutation. Use preset seed if provided.
    seed <- seeds[index]
    if(seed > 0)
      set.seed(seed[1])
    perms <- sample(seq(n.ind), n.ind, replace = FALSE)
    cross$pheno <- cross$pheno[perms,]
  }

  ## Cycle through all the phenotypes in sets of size size.set. Keeps R object smaller.
  for(pheno.set in seq(num.sets)) {
    cat(pheno.set, "\n")
    traits <- all.traits[trait.nums[trait.nums <= n.traits]]
    
    attach(file.path(dirpath, trait.file))
    perm.cross <- add.phenos(cross, get(trait.matrix)[, trait.names[traits], drop = FALSE],
                         index = trait.index)
    detach()

    pheno.col <- find.pheno(perm.cross, trait.names[traits])
    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar)

    per.scan.hl <- pull.highlods(per.scan, lod = lod.min, drop.lod = drop.lod)

    save(per.scan.hl, perms,
         file=file.path(dirpath, paste("per.scan",pheno.set, index,"RData",sep=".")))

    trait.nums <- trait.nums + size.set
  }

  chr.pos <- per.scan[,1:2]

  filenames <- list.files(dirpath, paste("per.scan", "[0-9][0-9]*", index, "RData",
                                         sep = "."))

  scan.hl <- cat.scanone(".", filenames, chr.pos)
  
  ## Has some problem with missing values for (any(diff(pos) < 0)) in runningmean.
  ## problem is in calc of maxlod.hl[[k]]: NA should be 0
  lod.sums <- lod.quantile.permutation(scan.hl,Nmax,lod.thrs,window,chr.pos,n.traits)
  
  save(scan.hl, lod.sums, perm.set,
       file = paste("perm", ".", cross.index, "_", index, ".RData", sep = ""))
}
#############################################################################################
big.phase3 <- function(dirpath = ".", index, ..., verbose = FALSE)
{
  ## Phase 3. Merge back together.

  load("Phase1.RData")

  n.thrs <- length(lod.thrs)
  
  max.lod.quant <- matrix(NA, n.perm, Nmax)
  max.N <- max.N.window <- matrix(NA, n.perm, n.thrs)
  
  for(i.perm in seq(n.perm)) {
    attach(file.path(dirpath, "perm", ".", cross.index, "_", i.perm, ".RData", sep = ""))
    n.quant <- length(lod.sums$max.lod.quant)
    max.lod.quant[i.perm, seq(n.quant)] <- lod.sums$max.lod.quant
    max.N[i.perm,] <- lod.sums$max.N$max.N
    max.N.window[i.perm,] <- lod.sums$max.N$max.N.win
  }
  save(max.lod.quant, max.N, max.N.window, file = "Phase3.RData", compress = TRUE)
}
