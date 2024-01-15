## Copyright 2016, 2017, 2018, 2022 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.


## Be VERY careful with what is done when frequencies are 0.  CBN will
## place events with frequency 0 without removing them.  Note which ones
## are removed in the code and then beware of what is reconstructed, as we
## will need to remove nodes but account for that later. The poset is
## always numbered 1 to number of nodes passed.  Done when reading
## posetToGraph, in the names argument.

## Note also that starting from an OT or similar is done in their 2009
## paper: see page 2812-2813, in the Jiang example.

## But then, in a sense it is being smart to place nodes with zero, or
## close to zero, freq. at the bottom of graph.

## Some functions limit columns/genes to 10. 12 was what Farahani and Lagergren
## did in theri comparison with H-CBN, and as Gerstung et al do.


## Testing CBN programs available
.._OncoSimul_test.ctcbn <- system("ct-cbn -h", ignore.stdout = TRUE)
.._OncoSimul_test.hcbn <- system("h-cbn -h", ignore.stdout = TRUE)
if(.._OncoSimul_test.ctcbn || .._OncoSimul_test.hcbn) {
    warning(paste(
        "\n\n",
        "\n******************************************************",
        "\n******************************************************\n",
        ## "          OncoSimulT installation warning:\n",
        "\n  WARNING \n",
        "The external programs h-cbn and ct-cbn were not found.",
        "You will not be able to use CBN.",
        "You can download them from http://www.bsse.ethz.ch/cbg/software/ct-cbn.html",
        "See further instructions in the README.",
        "\n******************************************************",
        "\n******************************************************",
        "\n\n"
        )
            )
}
rm(.._OncoSimul_test.ctcbn)
rm(.._OncoSimul_test.hcbn)


## be very careful using omp_threads > 1. I think OMP can actually
## slow things down.
cbn_proc <- function(x, addname, init.poset = "OT", nboot = 0,
                     verbose = FALSE, omp_threads = 1, silent = TRUE,
                     parall = FALSE, mc.cores = detectCores(),
                     rmfile = TRUE) {
    if(is.null(colnames(x))) stop("No colnames for input matrix")
    if(nboot == 0) {
        edges <- run.cbn(x, addname = addname, init.poset = init.poset,
                         omp_threads = omp_threads, silent = silent, rmfile = rmfile)
        edges <- cbind(edges, CBN_edgeBootFreq = NA,
                       stringsAsFactors = FALSE)
    } else {
      ## FIXME: we seem to always run with parall, but nboot = 0
      ## and we do not use mclapply.
      ## We should remove the !parall part of the if, or simplify this.
      ## deleting the !parall seems best, as the code in use now is more
      ## general and easy to extend
      if(!parall) { 
        ## Returns CBN and bootstrap freqs of edges
        edges <- run.cbn(x, addname = addname, init.poset = init.poset,
                         omp_threads = omp_threads, silent = silent, rmfile = rmfile)
        CBN_edgeBootFreq <- rep(0, nrow(edges))
        names(CBN_edgeBootFreq) <- edges[, "edge"]
        nn <- names(CBN_edgeBootFreq)
        for(i in seq.int(nboot)) {
          if(verbose)
            message("\n  .... doing bootstrap ", i, "\n")
          ind <- sample(nrow(x), nrow(x), replace = TRUE)
          bx <- x[ind, , drop = FALSE]
          addnameb <- paste0(addname, "b", i, paste(sample(letters, 3), collapse=""))
          bootedges <- run.cbn(bx, addname = addnameb, init.poset = init.poset,
                               omp_threads = omp_threads, silent = silent,
                               rmfile = rmfile)$edge
          posadd <- na.omit(match(bootedges, nn))
          if(length(posadd))
            CBN_edgeBootFreq[posadd] <- CBN_edgeBootFreq[posadd] + 1
        }
      } else {
        ## We parallelize over all calls, including the first,
        ## non-boottstrapped one
        ll <- c(0, seq.int(nboot))
        paralfunc <- function(index, odata = x,
                              addname_ = addname,
                              init.poset_ = init.poset, omp_threads_ = omp_threads,
                              silent_ = silent) {
          dd <- boot_data_index(odata, index)
          
          ## ## Save dd before cbn execution for debugging
          ## save(dd,file=paste0("dds/other/dd_",strsplit(ff[dbgfile],"\\.rds"),"_seed",toString(rseed),"_index",toString(index),".Rdata"))
          ## tmp <- run.cbn(dd, addname = "debug", init.poset = init.poset_,
          ##                omp_threads = omp_threads_, silent = TRUE, rmfile = TRUE)
          
          addnameb <- paste0(addname_, "b", index, paste(sample(letters, 3),
                                                         collapse=""))
          tmp <- run.cbn(dd, addname = addnameb, init.poset = init.poset_,
                         omp_threads = omp_threads_, silent = silent_,
                         rmfile = rmfile)
          
          if(index > 0) {
            return(tmp$edge)
          } else {
            return(tmp)
          }
        }
        ## Do not use mclapply at all. No parallel to allow multiple indep. runs
        alledges <- lapply(ll, function(z) paralfunc(z))
        edges <- alledges[[1]]
        CBN_edgeBootFreq <- rep(0, nrow(edges))
        names(CBN_edgeBootFreq) <- edges[, "edge"]
        nn <- names(CBN_edgeBootFreq)
        tedges <- table(unlist(alledges[-1])) ## gives all.
        nm_in_boot_table <- intersect(nn, names(tedges))
        if(length(nm_in_boot_table))
          CBN_edgeBootFreq[nm_in_boot_table] <- tedges[nm_in_boot_table]
      }
      edges <- cbind(edges, CBN_edgeBootFreq = CBN_edgeBootFreq/nboot,
                     stringsAsFactors = FALSE)
    }
    return(list(edges = edges,
                nboot = nboot,
                init.poset = init.poset))
}

run.cbn <- function(x,
                    init.poset = "linear", ## could be OT?
                    dirname = NULL,
                    addname = NULL, 
                    temp = 1,
                    ## as in Hosseini et al., p. i391
                    ## but for ncol(x) <= 4 200, not 100
                    steps = ifelse(ncol(x) <= 4, 200,
                            ifelse(ncol(x) == 5, 1000, 10000)),
                    ## max(200, ncol(x)^2), 
                    silent = TRUE,
                    eparam = 0.05, ## 0.05
                    rmfile = TRUE,
                    omp_threads = 1,
                    custom.poset = NULL
                    ) {

    if(is.null(steps))
        steps <- ncol(x)^2 ## Their default
    
    if(is.null(dirname)) {
        dirname <- tempfile()
        dirname0 <- NULL
        if(!is.null(addname)) {
            dirname0 <- dirname
            dirname <- paste0(dirname, "/",
                              "_cbn_", init.poset, "_",
                              addname)
        }
        if(!silent)
            message(paste("\n Using dir", dirname))
        if(dir.exists(dirname)) {
            stop("dirname ", dirname, "exists")
        }
        dir.create(dirname, recursive = TRUE)
    }

    
    ## their defaults are temperature = 1,
    ## and number of steps = number of genes ^ 2
    zzz <- call.external.cbn(x, createdir = FALSE,
                             dirname = dirname, eparam = eparam,
                             temp = temp, steps = steps,
                             silent = silent, init.poset = init.poset,
                             omp_threads = omp_threads, custom.poset = custom.poset)
    cnames <- colnames(x)
    poset <- read.poset(dirname, ncol(x))
    rm(zzz)
    ## Actually, leave the root in there. But when getting the transitive
    ## closure, do not use.  Why leave root? Because oncotree and treemix
    ## seem to use it. And because it is cleaner and clearer.

    gr <- evam_posetToGraph(poset, names = c("Root", cnames),
                                    addroot = TRUE, type = "adjmat",
                                    strictAdjMat = FALSE)
    gr <- sortAdjMat(gr)
    evam_checkProperMinimalAdjMat(gr, orderedNames = FALSE)
    
    lambdas <- read.lambda(dirname)
    rerun_lambda <- rerun_final_lambda(dirname, omp_threads = omp_threads)
    lambdas <- cbind(lambdas, rerun_lambda = rerun_lambda)
    
    rownames(lambdas) <- cnames
    
    if(rmfile) {
        files.created <- paste(dirname, c(".pat", ".prf", ".log", ".poset",
                                          ".lambda"), sep = "")
        file.remove(files.created)
        try(file.remove(paste0(dirname, ".lambda_i"))) ## for the hack that produces lambda_i
        try(file.remove(paste0(dirname, ".logliknew"))) ## for the 2nd hack that produces logliknew
        ## to get the lambda, iterative, file
        file.remove(paste(dirname, "/00000.poset", sep = ""))
        file.remove(dirname)
        if(!is.null(dirname0))
            file.remove(dirname0)
    }

    ## We simply want the edges list, as simpler for bootstrap later
    present <- which(gr == 1, arr.ind = TRUE)
    From <- rownames(gr)[present[, 1]]
    To <- colnames(gr)[present[, 2]]
    Edge <- paste(From, To, sep = " -> ")
    init_lambda <- lambdas[To, "init_lambda"]
    final_lambda <- lambdas[To, "final_lambda"]
    rerun_lambda <- lambdas[To, "rerun_lambda"]
    
    return(
        data.frame(From = From,
                   To = To,
                   edge = Edge,
                   init_lambda = init_lambda,
                   final_lambda = final_lambda,
                   rerun_lambda = rerun_lambda,
                   stringsAsFactors = FALSE))
}


## return as vector, say vout, which is a paste(From, " -> ", To).
## counts is rep(0, length(vout)) with names <- vout
## then do, in boot, where vboot is bootstrapped
##  count[ vboot %in% names(count) ] <- count[ vboot %in% names(count) ] + 1

## a brutish thing would be to accumulate the character vector and then do
## a table and subset.



sortAdjMat <- function(am) {
    ## If column names, except Root, are integers, sort as integers. O.w.,
    ## general lexicog. sort.
    cn <- colnames(am)
    rootpos <- grep("^Root$", cn) 
    if(length(rootpos) != 1)
        stop("No root in adj mat, or multiple Roots")
    cn0 <- colnames(am)[-rootpos]
    namesInts <- type.convert(cn0, as.is = TRUE)
    if(is.integer(namesInts)) {
        cn <- c("Root", evam_string_sort(namesInts))
    } else {
        cn <- c("Root", evam_string_sort(cn0))
    }
    return(am[cn, cn])
}



call.external.cbn <- function(data,
                              createdir = FALSE,
                              dirname = "testcbn", eparam = 0.05,
                              temp = 1, steps = 200, silent = TRUE,
                              init.poset = "linear", omp_threads = 1,
                              custom.poset = NULL) {
    ## I assume h-cbn and ct-cbn are available
    data2 <- cbind(1, data)
    write(c(nrow(data2), ncol(data2)),
          file = paste(dirname, ".pat", sep = ""),
          sep = " ")
    write(t(data2), file = paste(dirname, ".pat", sep = ""),
          ncolumns = ncol(data2),
          append = TRUE, sep = " ")

    ## Hangs.
    ## Reported: https://github.com/Rdatatable/data.table/issues/1727
    ## data.table::fwrite(data.table(data2),
    ##                    file = paste(dirname, ".pat", sep = ""),
    ##                    append = TRUE, sep = " ",
    ##                    col.names = FALSE)
    write(c("null", colnames(data)),
          file = paste(dirname, ".prf", sep = ""),
          sep = " ")
    if (is.null(omp_threads)) {
        ## setting this to detectCores() can lead to huge computing times
        ## if nothing passed, 1, just in case
        OMPthreads <- 1 
    } else{
        OMPthreads <- omp_threads
    }
    message("\n Exporting OMP_threads from call.external.cbn. OMP_NUM_THREADS = ", OMPthreads, "\n")
    ompt <- paste("export OMP_NUM_THREADS=", OMPthreads, "; ", sep = "")

    if(createdir)
        dir.create(dirname)
    if(!dir.exists(dirname))
        stop(" directory ", dirname, " does not exist")
    if(file.exists(paste0(dirname, "/00000.poset")))
        stop("A previous run in this directory?")
    if(file.exists(paste0(dirname, ".lambda")))
        stop("A previous run in this directory?")
    if(file.exists(paste0(dirname, ".poset")))
        stop("A previous run in this directory?")
    
    if(init.poset == "linear") {
        write.linear.poset(data, dirname)
    } else if (init.poset == "OT") {
        ##  will OT always return all of the nodes in the poset? Not
        ##  necessarily, but it is not a problem. CBN will start from the
        ##  poset, but explore around.
        ot1 <- try(run.oncotree(data)) ## , type.out = "adjmat"))
        if(inherits(ot1, "try-error")) {
            warning("\n   Error using OT for init poset; using linear.poset\n")
            write.linear.poset(data, dirname)
        } else {
            ## Remember that the posets always number starting from 1. But OT
            ## respects the names in the data. So we need to map those names
            ## to a sequence 1:number of columns in the data which is what the
            ## silly poset format understands.
            
            newnames <- c("Root", match(colnames(ot1)[-1], colnames(data)))
            colnames(ot1) <- rownames(ot1) <- newnames
            write.poset(evam_OTtoPoset(ot1),
                        ncol(data), dirname)
        }
    } ## The next are probably such bad ideas in almost
    ## any conceivable case that I comment the code

    ## else if (init.poset == "custom") {
    ##     write.poset(custom.poset, ncol(data), dirname)
    ## } else { ## Use ct-cbn to search and create starting poset;
    ##     ## possibly eternal. NOT RECOMMENDED
    ##     warning("Not using an initial poset can take VERY long")
    ##     writeLines(as.character(c(ncol(data), 0)),
    ##                con = paste(dirname, ".poset", sep = ""))
    ##     ## First create the lambda file
    ##     zzz <- system(paste(ompt , paste("h-cbn -f",  dirname, "-w")),
    ##                   ignore.stdout = silent)
    ##     if(!silent) cat("\n\n")
    ##     ## this call requires a lambda file
    ##     zzz <- system(paste(ompt,
    ##                         paste("h-cbn -f",  dirname,
    ##                               "-e", format(eparam, scientific = FALSE),
    ##                               "-w -m")), ignore.stdout = silent)
    ##     rm(zzz)
    ##     if(!silent) cat("\n\n")
    ## }
    
  ## Remove option -m, the printing of most likely path as
  ##    - we do not use it now
  ##    - it can lead to strange problems getting millions of ceros printed out
  
    zzz <- system(paste(ompt,
                        paste("h-cbn -f",  dirname, "-s", 
                              "-T", format(temp, scientific = FALSE),
                              "-N", format(round(steps), scientific = FALSE),
                              "-w")), ignore.stdout = silent)
    rm(zzz)
  if (!silent) cat("\n\n")
  ## the final poset in dirname/00000.poset
}

## For now, using my hack to get both initial and final
## lambda
read.lambda <- function(dirname, verbose = FALSE) {
    fname1 <- paste(dirname, ".lambda", sep = "") ## the first one
    fname2 <- paste(dirname, ".lambda_i", sep = "")
    init_lambda <- read.lambda.file(fname1, verbose = verbose)
    final_lambda <- read.lambda.file(fname2, verbose = verbose)
    return(cbind(init_lambda, final_lambda)[-1, ])
}

read.lambda.file <- function(fname, verbose = FALSE) {
    tmp <- try(scan(fname,
                    quiet = !verbose,
                    what = numeric()))
    if(inherits(tmp, "try-error")) stop("error in read.lambda.file",
                                        fname, "\n")
    return(tmp)
}

## re-run h-cbn to get lambda
rerun_final_lambda <- function(dirname, omp_threads, verbose = FALSE) {
    ## go to the dir
    ## cp dirname/00000.poset as dirname.poset
    ## call h-cbn -f dirname -w
    ## read the lambda file
    ## We are not setting OMP_NUM_THREADS here
    ## I guess we take the one in effect from the shell, it there is one
    ## In particular, the one set by former call to call.external.cbn
    ## is probably ignored?
    ## Now, when I launch many CPMs that use mclapply, I set
    ## OMP_NUM_THREADS = 1 in the bash script that calls R


    if (is.null(omp_threads)) {
        OMPthreads <- 1
    } else{
        OMPthreads <- omp_threads
    }
    message("\n Exporting OMP_threads from rerun_final_lambda. ",
            "OMP_NUM_THREADS = ", OMPthreads, "\n")
    ompt <- paste("export OMP_NUM_THREADS=", OMPthreads, "; ", sep = "")

    
    zzz <- system(paste0("cp ", dirname, "/00000.poset ",
                         dirname, ".poset"))
    zzz2 <- system(paste(ompt,
                         paste0("h-cbn -f ", dirname, " -w")),
                   ignore.stdout = !verbose)
    relambda <- read.lambda.file(paste(dirname, ".lambda", sep = ""),
                                 verbose = verbose)
    rm(zzz)
    rm(zzz2)
    return(relambda[-1])
}


read.poset.file <- function(dirname, maxn, verbose = FALSE) {
    ## Read a poset as generated by h-cbn
    tmp <- try(scan(dirname,
                quiet = !verbose,
                what = integer()))
    tmp0 <- scan(dirname,
                 quiet = !verbose)
    stopifnot(isTRUE(all.equal(tmp, tmp0)))
    if(inherits(tmp, "try-error")) stop("error in read.poset.file")
    tmp <- tmp[-(length(tmp))]
    nn <- tmp[1]
    tmp <- tmp[-1]
    tmp <- matrix(tmp, ncol = 2, byrow = TRUE)
    ## nodes with no ancestors or descendants
    missing.nodes <- setdiff(1:nn, unique(tmp))
    if(length(missing.nodes)) {
        if(verbose)
            message("Reading a poset with missing nodes") ## this is OK if a node not placed in the graph
        if(maxn != nn)
            stop("maxn != nn and missing nodes. Probably should not happen") ## FIXME: stop here?
        mnm <- matrix(c(missing.nodes, rep(NA, length(missing.nodes))),
                      ncol = 2) ## this is perfectly OK, and
                                ## poset.to.graph will do what it should
        tmp <- rbind(tmp, mnm)  ## ditto.
    } else {
        if(maxn != nn)  ## FIXME: I think I should stop in this case,
                        ## unless the user says o.w. But poset.to.graph will break anyway.
            stop("No missing nodes but maxn != nn. Probably should not happen")
    }
    return(tmp)
}

read.poset <- function(dirname, maxn, verbose = FALSE) {
    dirname <- paste(dirname, "/00000.poset", sep = "")
    read.poset.file(dirname, maxn = maxn, verbose = verbose)
}


linear.poset <- function(x) {
  ## a direct translation of linear_poset and write_poset
  ## in cbn.py by
  ## Niko Beerenwinkel and Moritz Gerstung

    ## nr is never used. The poset matrix is a square matrix
    ## nr <- nrow(x)
    ## See their code
    ## n, m = pat.shape
    ## poset = np.zeros([m,m])
  nc <- ncol(x)
  sorted <- order(colMeans(x), decreasing = TRUE)
  poset <- matrix(0, ncol = nc, nrow = nc)
  s <- sorted[1]
  for (t in sorted[2:nc]) {
    poset[s, t] <- 1
    s <- t
  }
 
  ## now, translate write_poset.
  ## posetw <- matrix(0, ncol = 2, nrow = nr)
  ## for (i in 1:nc) {
  ##   for(j in 1:nc) {
  ##     if(poset[i, j])
  ##       posetw[i, ] <- c(i, j)
  ##   }
  ## }
  
  ## do the R way
  posetw <- which(poset == 1, arr.ind = TRUE)
  posetw <- posetw[order(posetw[, 1]), ]  ## place drop = FALSE?
  return(posetw)
}


write.linear.poset <- function(x, dirname) {
    poset <- linear.poset(x)
    write.poset(poset, ncol(x), dirname)
}

write.poset <- function(poset, ncoldata, dirname) {
  fn <- paste(dirname, ".poset", sep = "")
  write(ncoldata, file = fn)
  write(t(poset), file = fn, append = TRUE,
        sep = " ", ncolumns = 2)
  write(0, file = fn, append = TRUE)
}


run.oncotree <- function(x, 
                         error.fun = function(x, y) { sum((x - y)^2)},
                         recover.failure = TRUE) {

    ## yes, ugly, but I do not want to specify it
    onco.fit <- try(oncotree.fit(x, error.fun = error.fun))
    if(inherits(onco.fit, "try-error") && recover.failure) {
            onco.fit <- try(oncotree.fit(x, error.fun = NULL))
    }
    if(inherits(onco.fit, "try-error")) {
        message("\n oncotree.fit failed")
        ## no, I do not return anything here. I want it to fail, and be
        ## captured with a try if needed.
    }
    
    gdf <- graph.data.frame(data.frame(parents = onco.fit$parent$parent[-1],
                                       children = onco.fit$parent$child[-1]
                                       ),
                            directed = TRUE,
                            vertices = NULL)
    ## This would be faster if getting a graphNEL, but unsafer so let's be
    ## very strict

    am <- get.adjacency(gdf, sparse = FALSE)
    am1 <- am
    storage.mode(am1) <- "integer"
    stopifnot(isTRUE(all.equal(am, am1)))
    sam1 <- sortAdjMat(am1)
    evam_checkProperOTAdjMat(sam1)
    return(sam1)
    ## if(type.out == "adjmat")
    ##     return(sam1)
    ## else if (type.out == "graphNEL")
    ##     return(as(sam1, "graphNEL"))
}
## library(codetools)
## checkUsageEnv(env = .GlobalEnv)
