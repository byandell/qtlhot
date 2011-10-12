big.phase0 <- function(dirpath, cross, trait.file, trait.matrix,
                       droptrait.names = NULL,
                       keeptrait.names = NULL,
                       lod.thrs,
                       sex = "Sex", trait.index, batch.effect = NULL, size.set = 250,
                       verbose = TRUE)
{
  ## Set up cross object and trait.data objects.

  ## Cross object with sex, trait.index and batch.effect as only phenotypes.
  traits <- match(c(sex,trait.index, batch.effect), names(cross$pheno), nomatch = 0)
  cross$pheno <- cross$pheno[, traits, drop = FALSE]
  
  ## Subset to individuals with batch if used.
  if(!is.null(batch.effect))
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))

  ## Save cross object and other pertinent information.
  if(verbose)
    cat("Saving cross object in cross.RData\n")
  save(cross, lod.thrs, file = "cross.RData", compress = TRUE)

  ## Get trait values for selected traits.
  attach(file.path(dirpath, trait.file))
  trait.names <- dimnames(get(trait.matrix))[[2]]
  all.traits <- seq(trait.names)
  if(!is.null(droptrait.names))
    all.traits[match(droptrait.names, trait.names, nomatch = 0)] <- 0
  if(!is.null(keeptrait.names))
    all.traits[-match(keeptrait.names, trait.names, nomatch = 0)] <- 0
  all.traits <- all.traits[all.traits > 0]
  n.traits <- length(all.traits)

  num.sets <- ceiling(n.traits / size.set)
  trait.nums <- seq(size.set)

  ## Save trait names, including those dropped, and size of sets.
  tmp <- paste("Trait", 0, "RData", sep = ".")
  if(verbose)
    cat("Saving trait metadata in", tmp, "\n")
  save(trait.file, trait.matrix, droptrait.names, size.set,
       trait.names, all.traits, num.sets,
       file = tmp, compress = TRUE)

  ## Cycle through all the phenotypes in sets of size size.set. Keeps R object smaller.
  for(pheno.set in seq(num.sets)) {
    traits <- all.traits[trait.nums[trait.nums <= n.traits]]
    
    trait.data <- get(trait.matrix)[, trait.names[traits], drop = FALSE]

    tmp <- paste("Trait", pheno.set, "RData", sep = ".")
    if(verbose)
      cat("Saving", length(traits), "trait data as set", pheno.set, "in", tmp, "\n")
    save(trait.data, pheno.set, file = tmp, compress = TRUE)

    ## Get next size.set traits.
    trait.nums <- trait.nums + size.set
  }
  detach()
}
#############################################################################################
big.phase1 <- function(dirpath = ".", cross.index = 0, params.file,
                       cross, lod.thrs, Nmax = 2000, n.perm = 1, n.split = 1,
                       batch.effect = NULL,
                       drop.lod = 1.5, lod.min = min(lod.thrs),
                       window = 5, seed = 0, big = TRUE, ...,
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

  ## Load trait.
  load(file.path(dirpath, "Trait.0.RData"))

  ## Random numbers.
  if(n.perm > 1) {
    if(seed[1] > 0)
      set.seed(seed[1])

    n.ind <- nind(cross)
    seeds <- sample(c(98765:987654321), n.perm, replace = FALSE)
  }
  else
    seeds <- 0

  ## Big assumes n.split is n.perm and does one perm each run.
  n.split <- n.perm
  
  save(cross, n.perm, seeds, Nmax, drop.lod, lod.thrs, lod.min, batch.effect,
       trait.index, trait.names, all.traits, size.set, ## From Trait.0.RData
       window, cross.index, n.split, big,
       file = "Phase1.RData", compress = TRUE)
  invisible()
}
#############################################################################################
big.phase2 <- function(dirpath = ".", index, ..., remove.files = TRUE, verbose = FALSE)
{
  ## Phase 2.
  ## Loop on sets of phenotypes, adding only what is needed.
  ## Loop on permutations internal to scanone.permutations.
  
  load(file.path(dirpath, "Phase1.RData"))

  if(verbose)
    cat("compute covariates\n")
  covars <- sexbatch.covar(cross, batch.effect, verbose = TRUE)

  n.ind <- nind(cross)
  n.phe <- nphe(cross)
  n.traits <- length(all.traits)

  if(n.perm > 1) {
    ## Random permutation. Use preset seed if provided.
    seed <- seeds[index]
    if(seed > 0)
      set.seed(seed[1])
    if(verbose)
      cat("sample permutation", seed[1], "\n")
    perms <- sample(seq(n.ind), n.ind, replace = FALSE)
    cross$pheno <- cross$pheno[perms,]
  }

  ## Cycle through all the phenotypes in sets of size size.set. Keeps R object smaller.
  ## Assume large trait matrix has been broken into Trait.i.RData, each with trait.data.

  filenames <- list.files(dirpath, "Trait.[1-9][0-9]*.RData")
  if(!length(filenames))
    parallel.error(5, 5, index)
  
  for(pheno.set in seq(length(filenames))) {
    if(verbose)
      cat(pheno.set, "\n")

    ## This is not working correctly. Either Trait.n.RData are all the same or ??.
    attach(file.path(dirpath, filenames[pheno.set]))
    perm.cross <- add.phenos(cross, trait.data, index = trait.index)
    pheno.col <- find.pheno(perm.cross, dimnames(trait.data)[[2]])
    detach()

    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar)

    per.scan.hl <- pull.highlods(per.scan, lod = lod.min, drop.lod = drop.lod)

    save(per.scan.hl, perms,
         file=file.path(dirpath, paste("per.scan",pheno.set, index,"RData",sep=".")))
  }

  chr.pos <- per.scan[,1:2]

  filenames <- list.files(dirpath, paste("per.scan", "[0-9][0-9]*", index, "RData",
                                         sep = "."))

  scan.hl <- cat.scanone(dirpath, filenames, chr.pos)

  if(remove.files)
    file.remove(dirpath, filenames)
  
  ## Has some problem with missing values for (any(diff(pos) < 0)) in runningmean.
  ## problem is in calc of maxlod.hl[[k]]: NA should be 0
  lod.sums <- lod.quantile.permutation(scan.hl,Nmax,lod.thrs,window,chr.pos,n.traits)
  
  save(scan.hl, lod.sums, cross.index, index,
       file = paste("perm", ".", cross.index, "_", index, ".RData", sep = ""))
}
#############################################################################################
big.phase3 <- function(dirpath = ".", index, cross.index, ..., verbose = FALSE)
{
  ## Phase 3. Merge back together.

  load(file.path(dirpath, "Phase1.RData"))

  n.thrs <- length(lod.thrs)
  
  max.lod.quant <- matrix(NA, n.perm, Nmax)
  max.N <- max.N.window <- matrix(NA, n.perm, n.thrs)
  
  for(i.perm in seq(n.perm)) {
    attach(file.path(dirpath, paste("perm", ".", cross.index, "_", i.perm, ".RData", sep = "")))
    n.quant <- length(lod.sums$max.lod.quant)
    max.lod.quant[i.perm, seq(n.quant)] <- lod.sums$max.lod.quant
    max.N[i.perm,] <- lod.sums$max.N$max.N
    max.N.window[i.perm,] <- lod.sums$max.N$max.N.win
    detach()
  }
  outfile <- paste("Phase3", ifelse(cross.index > 0, cross.index, ""), ".RData", sep = "")

  save(max.lod.quant, max.N, max.N.window, file = outfile, compress = TRUE)
}
