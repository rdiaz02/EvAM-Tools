# Copyright 2022 Ramon Diaz-Uriarte, Pablo Herrera Nieto.

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## ## timings
## dates_for_timing <- function(x) {
##     cat("\n At  ;", x, ";", date(), "\n")
## }






## Sample an indivial from a transition rate matrix}
## \item{trm}{transition rate matrix}
## \item{T_sampling}{Time at which sampling happens.}
## \item{ngenots}{Number of genotypes}
## \item{genot_names}{String array with genotype names}
## \value{
## sampled genotype, trajectory, and accumulated time
## }
## \description{
## This is a transition *rate* matrix, *not* a matrix
## of transition probabilities between genotypes. For example, this
## cannot be used with OT.
## }

indiv_sample_from_trm <- function(trm, T_sampling, ngenots = NULL,
                            genot_names = NULL) {
    if(is.null(ngenots)) ngenots <- ncol(trm)
    if(is.null(genot_names)) genot_names <- colnames(trm)
    row <- 1
    t_accum <- 0
    genotype <- "WT" 
    trajectory <- "WT"

    while(TRUE) {
        qii <- sum(trm[row, ]) ## Or qs in Gotovos et al terminology
        ## Special case of genotype that does not transition to anything
        if(qii == 0.0) break 
        t_transition <- rexp(n = 1, rate = qii)
        t_accum <- t_accum + t_transition
        if(t_accum >= T_sampling ) break
        ## We transition. To what genotype?
        pj <- trm[row, ]/qii
        j <- sample(x = seq_len(ngenots), size = 1, prob = pj)
        genotype <- genot_names[j]
        trajectory <- c(trajectory, genotype)
        ## For next iteration
        row <- j
    }
    return(list(genotype = genotype,
                trajectory = trajectory,
                t_accum = t_accum))
}

##
## @title Sample an individual from a transition rate matrix with precomputed
## matrix diagonal (and number of genotypes and genotype names).
##
## @description If we sample many individuals from the same process, it is worth
##     passing some common, precomputed values. This is what I call
##     "transition matrix standardized": Diagonal is passed separately, and
##     entries in the matrix are probabilities.
## 
## This is all for a transition *rate* matrix, *not* a matrix
## of transition probabilities between genotypes. For example, this
## cannot be used with OT.
##
## @param trmstd transition rate matrix "standardized",
## @param diag diagonal of transition rate matrix, time of sampling of a case/individual
## @param T_sampling Time at which sampling happens.
## @param ngenots Number of genotypes
## @param genot_names String array with genotype names
## 
## @return sampled genotype, trajectory, and accumulated time
indiv_sample_from_trm_pre <- function(trmstd,
                                      diag,
                                      T_sampling,
                                      ngenots,
                                      genot_names) {
    row <- 1
    t_accum <- 0
    genotype <- "WT" 
    trajectory <- "WT"
    
    while (TRUE) {
        ## qii <- diag[row] ## Or qs in Gotovos et al terminology
        ## Special case of genotype that does not transition to anything
        if (diag[row] == 0.0) break 
        t_transition <- rexp(n = 1, rate = diag[row])
        t_accum <- t_accum + t_transition
        if(t_accum >= T_sampling ) break
        ## We transition. To what genotype?
        j <- sample(x = seq_len(ngenots), size = 1, prob = trmstd[row, ])
        genotype <- genot_names[j]
        trajectory <- c(trajectory, genotype)
        ## For next iteration
        row <- j
    }
    return(list(genotype = genotype,
                trajectory = trajectory,
                t_accum = t_accum))
}

##
## @title Sample a population from a transition rate matrix.
## 
## @description Like indiv_sample_from_trm, but for multiple individuals.
## 
## @param trm transition rate matrix, number of samples or times of samples,
## @param n_samples Int with the number of samples to be computed
## @param T_sampling Time at wich each individual in sample. By default they 
## are randomly generated from an exponential distribution of rate 1.
## @param pre_compute whether or not to precompute entries of the trans rate matrix (for speed)
## @param cores number of cores (pass 1 is you do not want to parallelize)
## 
## @return List with precompunted sampling time of sampling, the actual time sampled
## observed for each sample (s), the complete trajectory of acquired mutations
## and the observed genotype

population_sample_from_trm <- function(trm, n_samples = 10,
                                       T_sampling = NULL,
                                       pre_compute = TRUE,
                                       cores = detectCores()) {
    if(is.null(T_sampling) && is.null(n_samples)) {
        stop("Pass either n_samples or T_sampling vector")
    }
    if(!is.null(T_sampling) && !is.null(n_samples)) {
        message("Ignoring n_samples as passing T_sampling")
    }
    if(is.null(T_sampling) && !is.null(n_samples)) {
        T_sampling <- rexp(n = n_samples, rate = 1)
    }

    ngenots <- ncol(trm)
    genot_names <- colnames(trm)
    if(pre_compute) {
        ## Like code in trans_rate_to_trans_mat
        sx <- rowSums(trm)
        ii <- which(sx > 0)
        for(i in ii) {
            trm[i, ] <- trm[i, ]/sx[i]
        }
        
        out <- mclapply(T_sampling,
               function(x)
                   indiv_sample_from_trm_pre(trmstd = trm,
                                             diag = sx,
                                             T_sampling = x,
                                             ngenots = ngenots,
                                             genot_names = genot_names),
               mc.cores = cores)  
        
    } else {
        out <- mclapply(T_sampling,
               function(x)
                   indiv_sample_from_trm(trm = trm,
                                         T_sampling = x,
                                         ngenots = ngenots,
                                         genot_names = genot_names),
               mc.cores = cores)
    }
    ## Structure output as Pablo's  simulate_population_2
    ## Otherwise, we could just exit from the above
    ## This will add time and increase RAM usage
    return(list(
        T_sampling = T_sampling
      , T_sum_events = unlist(lapply(out, function(x) x$t_accum))
      , trajectory = lapply(out, function(x) x$trajectory)
      , obs_events = unlist(lapply(out, function(x) x$genotype))
        ))
}

## # @noRd
## #' @title Process samples
## #' 
## #' @description Generate trajectories from data simulated from a given model.
## #' 
## #' @param sim list generated with population_sample_from_trm. Relevant
## #' fields are described below
## #' $T_sum_events time of events for the mutations of each gene
## #' $obs_events data.frame with mutated before the end of the sampling time
## #' @param n_genes number of genes observed
## #' @param gene_names List of gene names. Required.
## #' @param output type of output that we want
## #' 
## #' @return List with a list of trajectories (the order in which gene mutations
## #' are acquired), genotype frequencies and genotypes transition matrix (with
## #' counts of how many transitions between each genotype have been observed) 
process_samples <- function(sim, n_genes,
                            gene_names,
                            output = c("sampled_genotype_counts",
                                       "state_counts",
                                       "obs_genotype_transitions"),
                            cores = detectCores()) {

    #Checking input
    params <- c("trajectory", "obs_events")
    for (i in params){
        if (!(i %in% names(sim))) 
            stop(sprintf("%s is missing from your samples", i))
    }

    ## Checking output variables
    valid_output <- c("sampled_genotype_counts",
                      "state_counts",
                      "obs_genotype_transitions")

    if (any(!(output %in% valid_output ))) stop("Incorrect output specified")
    if (length(output) == 0) stop("No output specified")
    
    ## Set up
    retval <- list()
    ## n_states <- 2^n_genes
    ## sorted_genotypes <- generate_pD_sorted_genotypes(n_genes, gene_names)

    ## Calculate frequencies: genotype frequencies
    if ("sampled_genotype_counts" %in% output) {
        retval$sampled_genotype_counts <-
            sample_to_named_pD_ordered_out(sim$obs_events,
                                           n_genes, gene_names, "vector")
    }

    ## Calculate transitions

    ## But really, why do we want this? The expected number of transitions is
    ## already known. What does this give us? A way of testing, but that is
    ## all. To see the sampling variability?
    
    ## This assumes that genotypes, as given by
    ## sorted_genotypes, correspond to the genotypes as they exist
    ## in sim$trajectory (thus sim$obs_events).
    ## sorted_genotypes uses sorting
    ## of gene names. But the transition rate matrices might not unless
    ## they have been computed that way. They are now, though.

    ## This is the slowest part. Two implementations. 
    if ("obs_genotype_transitions" %in% output) { ## observed genotype transitions
        unlisted_trajectories <- unlist(sim$trajectory)
        ## ## Implementation 1
        ## tt <- sparse_transM_from_genotypes(unique(unlisted_trajectories))
        ## for (traj in sim$trajectory) {
        ##     steps <- length(traj) - 1
        ##     if (steps > 0) {
        ##         for (i in seq_len(steps)) {
        ##             tt[traj[i], traj[i + 1]] <- tt[traj[i], traj[i + 1]]  + 1
        ##         }
        ##     }
        ## }
        ########      Implementation 2. Can be much faster
        ##   Create a matrix of indexes first, then count cases
        ##   and assign the entries in the sparse matrix
        ##   Creating  matrix of indeces can use mclapply

        
        tt2 <- sparse_transM_from_genotypes(unique(unlisted_trajectories))
        tindex <-  seq_along(rownames(tt2))
        names(tindex) <- rownames(tt2)
        get_index_pairs <- function(v, thetindex = tindex) {
            lv <- length(v) - 1
            if (lv > 0) {
                w <- unname(thetindex[v])
                cbind(w[seq_len(lv)], w[2:(lv + 1)])
            }
        }
        ## Get matrix of indeces of transitions
        all_pairs2 <- mclapply(sim$trajectory,
                               get_index_pairs,
                               mc.cores = cores)
        all_pairs_stacked <- do.call(rbind, all_pairs2)


        ## Count the cases of each index. Sort. Then traverse and add.
        ## When change in state, assign value to sparse matrix.
        ##   It could possibly be made faster using RLE on ordered pairs.
        oo <- order(all_pairs_stacked[, 1], all_pairs_stacked[, 2])
        apso <- all_pairs_stacked[oo, ]

        last <- nrow(apso)
        sum <- 1
        for(k in seq_len(last - 1)) {
            if (all(apso[k, ] == apso[k + 1, ])) {
                sum <- sum + 1
            } else {
                tt2[apso[k, 1], apso[k, 2]] <- sum
                sum <- 1
            }
        }
        ## There is a single last one
        if (!(all(apso[last - 1, ] == apso[last, ])))
            tt2[apso[last, 1], apso[last, 2]] <- 1
        else ## no change between last two ones
            tt2[apso[last, 1], apso[last, 2]] <- sum
        
        ## stopifnot(identical(tt, tt2))
        retval$obs_genotype_transitions <- tt2
    }


    ## Calculate state_counts
    if ("state_counts" %in% output) { ## times each genotype visited
        if (!("obs_genotype_transitions" %in% output))
            unlisted_trajectories <- unlist(sim$trajectory)

        retval$state_counts <-
            sample_to_named_pD_ordered_out(unlisted_trajectories,
                                           n_genes, gene_names, "data.frame")

        ## FIXME: Paranoid check. Will remove later
        if ("obs_genotype_transitions" %in% output) {
            cstt2 <- colSums(tt2)[-1]
            ## No WT
            statecounts_vector <- retval$state_counts[-1, "Counts"]
            names(statecounts_vector) <- retval$state_counts[-1, "Genotype"]
            statecounts_vector <- statecounts_vector[statecounts_vector > 0]
            cstt2 <- cstt2[order(names(cstt2))]
            statecounts_vector <- statecounts_vector[order(names(statecounts_vector))]
            stopifnot(isTRUE(all(cstt2 == statecounts_vector)))
        }
    }
    return(retval)
}

## Create empty sparse matrix from vector of genotypes
## with rows and columns ordered with WT first, and
## then by increasing mutations and within number of mutations, sorted
sparse_transM_from_genotypes <- function(genots) {
    genots <- unique(genots)
    muts <- stringi::stri_count_fixed(genots, ",")
    names(muts) <- genots
    muts["WT"] <- -1
    genots <- genots[order(muts, genots)]
    return(sparseMatrix(i = NULL, j = NULL,
                        dims = c(length(genots), length(genots)),
                        x = 0,
                        dimnames = list(genots, genots)))
}

## Return vector of ranks (not order) so we can generate
## an "Index" column, or similar, so that ordering according to it
## gives us standard order.
## x: genotype names
standard_rank_genots_1 <- function(x) {
    muts <- stringi::stri_count_fixed(x, ",")
    xwt <- which(x == "WT")
    if (length(xwt)) muts[xwt] <- -1
    return(order(order(muts, x)))
}

## Like version _1, but for two sets of genotypes, such as "From" and "To"
standard_rank_genots_2 <- function(x, y) {
    mutsx <- stringi::stri_count_fixed(x, ",")
    xwt <- which(x == "WT")
    if (length(xwt)) mutsx[xwt] <- -1
    mutsy <- stringi::stri_count_fixed(y, ",")
    ywt <- which(y == "WT")
    if (length(ywt)) mutsy[ywt] <- -1
    return(order(order(mutsx, x, mutsy, y)))
}


sample_CPMs <- function(cpm_output
                      , N
                      , methods = c("OT", "OncoBN",
                                    "CBN", "MCCBN",
                                    "MHN", "HESBCN")
                      , output = c("sampled_genotype_counts")
                      , obs_noise = 0
                      , genotype_freqs_as_data = TRUE
                        ) {
    retval <- list()

    ## Anything that is passed as "cpm_output" must have
    ## something_predicted_genotype_freqs
    methods <- unique(methods)
    available_methods <- vapply(methods,
                                function(m) {
                                    outn <- paste0(m, "_predicted_genotype_freqs")
                                    return(!(is.null(cpm_output[[outn]])) &&
                                           !(is.na(cpm_output[outn])))
                                }, logical(1)
                                )
    available_methods <- names(which(available_methods))
    
    l_methods <- length(available_methods)
    if (l_methods < length(methods)) {
        warning("At least one method you asked to be sampled ",
                "did not have output.")
    }
    methods <- available_methods
    
    output <- unique(output)
    valid_output <- c("sampled_genotype_counts",
                      "obs_genotype_transitions",
                      "state_counts")
    not_valid_output <- which(!(output %in% valid_output))
    if (length(not_valid_output)) {
        warning("Output(s) ",
                paste(output[not_valid_output], sep = ", ", collapse = ", "),
                " not among the available output.",
                " Ignoring the invalid output.")
        output <- output[-not_valid_output]
    }
    if (length(output) == 0) stop("No valid output given.")
    if (any(c("state_counts", "obs_genotype_transitions") %in% output)) {
        message("For the requested output we will need to simulate ",
                "from the transition rate matrix.")
    }

    some_pred <- cpm_output[[paste0(methods[1], "_predicted_genotype_freqs")]]
    gene_names <- sort(setdiff(unique(unlist(strsplit(names(some_pred),
                                                      split = ", "))),
                          "WT"))
    n_genes <- length(gene_names)

    for (method in methods) {
        if (method %in% c("OT", "OncoBN")) {
            genots_pred <- cpm_output[[paste0(method, "_predicted_genotype_freqs")]]

            retval[[sprintf("%s_sampled_genotype_counts", method)]] <-
                genot_probs_2_pD_ordered_sample(genots_pred,
                                                n_genes, gene_names, N,
                                                "vector")
            retval[[sprintf("%s_obs_genotype_transitions", method)]] <- NULL
        } else {
            ## evam always returns the method_trans_rate_mat
            ## even if just with an NA
            if (any(c("state_counts", "obs_genotype_transitions") %in% output)) {
                ## Need to simulate from trm
                trm <- cpm_output[[sprintf("%s_trans_rate_mat", method)]]
                ## Do we have the necessary output to simulate?
                if ((length(trm) == 1) && is.na(trm)) {
                    retval[[sprintf("%s_obs_genotype_transitions", method)]] <- NULL

                    ogt <- cpm_output[[sprintf("%s_predicted_genotype_freqs",
                                               method)]]
                    if ((length(ogt) == 1) && is.na(ogt)) {
                        retval[[sprintf("%s_sampled_genotype_counts", method)]] <- NULL
                    } else {
                        ## Yes, we could sample. But this should never happen.
                        stop("No transition rate matrix in output ",
                             "but predicted_genotype_freqs")
                    }
                } else { ## transition rate matrix present
                    sims <- population_sample_from_trm(trm, n_samples = N)

                    psamples <-
                        process_samples(sims,
                                        n_genes,
                                        gene_names,
                                        output = output)
                    
                    for (reqout in valid_output) {
                        if (reqout %in% output)
                            retval[[paste0(method, "_", reqout)]] <-
                                psamples[[reqout]]
                        else
                            retval[[paste0(method, "_", reqout)]] <- NA
                    }
                }
            } else {
                ## Multinomial sampling from the predicted genotypes
                ## as for OT and OncoBN
                genots_pred <- cpm_output[[sprintf("%s_predicted_genotype_freqs",
                                                   method)]]
                if ((length(genots_pred) == 1) && is.na(genots_pred)) {
                    retval[[sprintf("%s_sampled_genotype_counts",
                                    method)]] <- NULL
                } else {
                    retval[[sprintf("%s_sampled_genotype_counts", method)]] <-
                        genot_probs_2_pD_ordered_sample(genots_pred,
                                                        n_genes,
                                                        gene_names,
                                                        out = "vector")
                }
            }
        }
        if ("sampled_genotype_counts" %in% output) {
            ## #######################################
            ##
            ## Noise and genotype frequencies as data
            ##
            ## #######################################

            ## obs_noise == 0  && !genotype_freqs_as_data : do nothing
            ## obs_noise == 0  && genotype_freqs_as_data  : return also the data
            ## obs_noise >  0  && !genotype_freqs_as_data :
            ##                           - generate data
            ##                           - add noise to data
            ##                           - data as freqs and overwrite the object
            ##                           - rm data
            
            ## obs_noise >  0  && genotype_freqs_as_data : as former,
            ##                                             returning data

            data_name <- paste0(method, "_sampled_genotype_counts_as_data")
            gf_name   <- paste0(method, "_sampled_genotype_counts")
            if ((obs_noise == 0) && (genotype_freqs_as_data)) {
                retval[[data_name]] <-
                    genotypeCounts_to_data(retval[[gf_name]], e = 0)
            } else if (obs_noise > 0) {
                data_noised <- genotypeCounts_to_data(retval[[gf_name]],
                                                      e = obs_noise)
                retval[[gf_name]] <- reorder_to_pD(data_to_counts(data_noised,
                                                                  out = "vector"))

                if (genotype_freqs_as_data) {
                   retval[[data_name]] <- data_noised
                }
            }
        }
    }
    return(retval)
}

## @title Count genotypes 
## 
## Take a sample (a vector), with genotypes as "A, B", etc
## and return a vector of frequencies (counts) in the exact same
## order as used by MHN
## A much faster implementation

## @param x vector of genotypes
## @param ngenes total number of genes
## @param gene_names List of gene names. If NULL, genes will be named alphabetically
## 
## @return counts of all genotypes in same order as used by MHN
##         gene_names is always sorted inside the function to ensure
##         results do not differ by gene_names order


## Why do you need to provide ngenes and/or gene_names?
##   Because you might have taken a small sample where not all
##   genes you expect are represented.
sample_to_pD_order <- function(x, ngenes, gene_names = NULL) {
    if (is.null(gene_names)) gene_names <- LETTERS[seq_len(ngenes)]
    ## Consistency checks
    stopifnot(ngenes == length(gene_names))
    gene_names_in_x <- setdiff(unique(unlist(strsplit(unique(x), ", "))), "WT")
    stopifnot(length(gene_names_in_x) <= ngenes)
    stopifnot(all(gene_names_in_x %in% gene_names))
    rm(gene_names_in_x)

    x <- as.data.frame(table(x), stringsAsFactors = FALSE)
    gene_names <- sort(gene_names)
    genot_int <- x[, 1]
    genot_int <- gsub("^WT$", "", genot_int, fixed = "FALSE")

    genot_int <-
        vapply(genot_int,
               function(z)
                   State.to.Int(
                       as.integer(gene_names %in%
                                  strsplit(z, ", ", fixed = TRUE)[[1]])),
               numeric(1))

    return(tabulate(rep(unname(genot_int), x[, 2]),
                   nbins = 2^ngenes))
}


## Given a named vector, where names are genotypes,
## return the same vector sorted in the same order as sample_to_pD
##  i.e., in the same order that MHN and generate_pD_sorted_genotypes
##  use.
## Assumption: genes are ONLY those present in the set of genotypes
##  no genes that always were absent.
##  Missing genotypes left as NA
reorder_to_pD <- function(x) {
    genots_n <- names(x)
    ## Ensure genotype names are canonical
    genots_n <- canonicalize_genotype_names(genots_n)
    names(x) <- genots_n

    ## Get gene names
    gene_n <- unique(stringi::stri_replace_all_fixed(
                                  unlist(
                                      stringi::stri_split_fixed(genots_n, ",")),
                                  " ", ""))
    gene_n <- setdiff(gene_n, "WT")
    sorted_genots <- generate_pD_sorted_genotypes(length(gene_n),
                                               gene_names = gene_n)

    if (!all(genots_n %in% sorted_genots))
        stop("At least one genotype name not in sorted_genots")

    ## Sort to original, return
    x2 <- x[sorted_genots]
    names(x2) <- sorted_genots
    return(x2)
}

## given a vector of genotype names
##  - sort gene names
##  - separate gene names by ", "
canonicalize_genotype_names <- function(x) {
    no_space <- stringi::stri_replace_all_regex(x, pattern = "[\\s]", "")

    pasted_sorted <- unlist(lapply(strsplit(no_space, split = ",", fixed = TRUE),
                  function(v) (paste(sort(v), collapse = ", "))))
    return(pasted_sorted)
}

## Genotypes in what for me is their "standard, sensible, order"
## By number of mutations, and within number of mutations, ordered
## as given by order.
genotypes_standard_order <- function(gene_names) {
    gene_names <- sort(gene_names)
    allgt <- allGenotypes_3(length(gene_names))$mutated
    gtn <- vapply(allgt, function(v) paste(gene_names[v], collapse = ", "),
                  "")
    gtn[1] <- "WT"
    return(gtn)
}

## Given a named vector, where names are genotypes,
## return the same vector sorted in the same order as genotypes_standard_order
## Like reorder_to_pD, but in what for me is much more sensible
reorder_to_standard_order <- function(x) {
    genots_n <- names(x)
    ## Ensure genotype names are canonical
    genots_n <- canonicalize_genotype_names(genots_n)
    names(x) <- genots_n

    ## Get gene names
    gene_n <- unique(stringi::stri_replace_all_fixed(
                                  unlist(
                                      stringi::stri_split_fixed(genots_n, ",")),
                                  " ", ""))
    gene_n <- setdiff(gene_n, "WT")
    sorted_genots <- genotypes_standard_order(gene_names = gene_n)

    if (!all(genots_n %in% sorted_genots))
        stop("At least one genotype name not in sorted_genots")

    ## Sort to original, return
    x2 <- x[sorted_genots]
    names(x2) <- sorted_genots

    return(x2)
}



## Given a matrix of 0/1 return the genotypes in canonicalized way
genot_matrix_2_vector <- function(x) {
    gn <- colnames(x)
    gt1 <- apply(x, 1, function(v) paste(sort(gn[v == 1]), collapse = ", "))
    gt1[gt1 == ""] <- "WT"
    return(gt1)
}


## Obtain probabilities of genotypes from transition rate matrix
## under sampling time distributed as exponential rate 1.
## 
## Using equation 4 (p. 243) in Schill et al., 2020, Bioinformatics, 36
##    "Modelling cancer progression using Mutual Hazard Networks"
##    and following their code (but using Jacobi from Rlinsolve).
##
##    Assumptions:
##     - x is a sparse matrix
##     - First column/row of x is WT
##     - For now, the initial distribution is 100% are WT
##          (could change, ensuring genotype order matches)

##  The final all.equal uses a tolerance larger than that of
##  the usual all.equal.

## Yes, this is much slower, like up to one order of magnitude,
## than Schill's Generate.pTh. Still, it takes generally < 0.04 seconds
## for 4 to 8 genes, and most of it is spent in the checks.
probs_from_trm <- function(x,
                           tolerance = 10 * sqrt(.Machine$double.eps),
                           all_genotypes = TRUE) {
    p0 <-  c(1, rep(0, nrow(x) - 1))

    ## Recall our trans. rate matrix are. rows: from; columns: to.
    ## But in Schill they are transposed.
    if (Matrix::nnzero(tril(x)))
        stop("Lower triangular not 0. Is this transposed?")
     Q <- t(x)
    diag(Q) <- -1 * colSums(Q)

    ## Equation 4 in Schill et al. Thus
    ## (I - Q) * p = p0
    I_Q <- Matrix::Diagonal(nrow(Q)) - Q

    ## Limited experiments showed that for nrow(Q) < 1024
    ## fastmatrix's Gauss-Seidel was a lot faster, and consistently so.
    ## For larger, I guess Rlinsolve's usage of sparse matrices
    ## gives an edge. Why is Jacobi faster? No idea; parallel updates?
    ## (and Gauss-Seidel, from Rlinsolve, was sometimes faster but sometimes
    ##  much slower)
    if (nrow(Q) >= 1024) {
        p2 <- Rlinsolve::lsolve.jacobi(A = I_Q, B = p0, adjsym = FALSE,
                                       reltol = 1e-5 * sqrt(.Machine$double.eps),
                                       weight = 1, verbose = FALSE)
        p <- as.vector(p2$x)
    } else {
        p4 <- fastmatrix::seidel(a = as.matrix(I_Q), b = p0,
                                 tol = 1e-5 * sqrt(.Machine$double.eps),
                                 maxiter = 1000,
                                 start = rep(0, nrow(Q)))
        p <- as.vector(p4)
    }

    names(p) <- colnames(x)
    if (!isTRUE(all.equal(sum(p), 1))) {
        warning("sum(p) - 1  = ", sum(p) - 1)
    }

    if (!all_genotypes) return(p)

    ## Get genes and from them number genotypes and identity genotypes
    gene_names <- sort(setdiff(unique(unlist(strsplit(colnames(x), split = ", "))),
                          "WT"))
    number_genes <- length(gene_names)
    num_genots <- 2^number_genes

    if (length(p) == num_genots) return(p)

    allGts <- genotypes_standard_order(gene_names)
    p_all <- rep(0.0, length = length(allGts))
    names(p_all) <- allGts
    p_all[names(p)] <- p
    return(p_all)
}


random_evam <- function(ngenes = NULL, gene_names = NULL,
                                 model = c("OT", "CBN", "HESBCN", "MHN",
                                           "OncoBN")
                               , graph_density = 0.35
                               , cbn_hesbcn_lambda_min = 1/3
                               , cbn_hesbcn_lambda_max = 3
                               , hesbcn_probs = c("AND" = 1/3,
                                                  "OR" = 1/3,
                                                  "XOR" = 1/3)
                               , ot_oncobn_weight_min = 0
                               , ot_oncobn_weight_max = 1
                               , ot_oncobn_epos = 0.1
                               , oncobn_model = "DBN"
                                 ) {
    stopifnot((graph_density <= 1) && (graph_density >= 0))
    if (!(xor(is.null(ngenes), is.null(gene_names))))
        stop("Give exactly one of ngenes XOR gene_names")

    if (length(model) > 1) {
        warning("Only one model should be specified. Using the first")
        model <- model[1]
    }

    if (model == "MCCBN") {
        message("Generating a random MCCBN model is the same ",
                "as generating a random CBN model. We'll do that.")
        model <- "CBN"
    }

    if (is.null(gene_names)) gene_names <- LETTERS[seq_len(ngenes)]
    if (is.null(ngenes)) ngenes <- length(gene_names)
    gene_names <- sort(gene_names)
    output <- list()
    if (model == "MHN") {
        mhn_sparsity <- 1 - graph_density
        ## Recall these thetas have theta_i,j: effect of j on i.
        thetas <- Random.Theta(n = ngenes, sparsity = mhn_sparsity)
        colnames(thetas) <- rownames(thetas) <- gene_names
        output <- MHN_from_thetas(thetas)
    } else if (model == "CBN") {
        poset <- mccbn::random_poset(ngenes, graph_density = graph_density)
        lambdas <- runif(ngenes, cbn_hesbcn_lambda_min, cbn_hesbcn_lambda_max)
        names(lambdas) <-  colnames(poset) <- rownames(poset) <- gene_names
        output <- CBN_from_poset_lambdas(poset, lambdas)
    } else if (model == "HESBCN") {
        hesbcn_probs <- hesbcn_probs/sum(hesbcn_probs)
        poset <- mccbn::random_poset(ngenes, graph_density = graph_density)
        lambdas <- runif(ngenes, cbn_hesbcn_lambda_min, cbn_hesbcn_lambda_max)
        names(lambdas) <-  colnames(poset) <- rownames(poset) <- gene_names
        stopifnot(identical(sort(names(hesbcn_probs)), c("AND", "OR", "XOR")))
        output <-
            HESBCN_from_poset_lambdas_relation_probs(poset,
                                                     lambdas,
                                                     hesbcn_probs)
    } else if (model == "OT") {
        poset <- OT_random_poset(ngenes,
                                 graph_density = graph_density)
        weights <- runif(ngenes, ot_oncobn_weight_min, ot_oncobn_weight_max)
        names(weights) <-  colnames(poset) <- rownames(poset) <- gene_names
        output <- OT_from_poset_weights_epos(poset, weights, ot_oncobn_epos)

    } else if (model == "OncoBN") {
        poset <- mccbn::random_poset(ngenes,
                                     graph_density = graph_density)
        thetas <- runif(ngenes, ot_oncobn_weight_min, ot_oncobn_weight_max)
        names(thetas) <-  colnames(poset) <- rownames(poset) <- gene_names
        output <- OncoBN_from_poset_thetas_epsilon_model(poset, thetas,
                                                       ot_oncobn_epos,
                                                       oncobn_model)
    }

    if (model %in% c("CBN", "MHN", "HESBCN")) {
        outname <- paste0(model, "_predicted_genotype_freqs")
        inname  <- paste0(model, "_trans_rate_mat")
        output[[outname]] <- probs_from_trm(output[[inname]])
    }
    return(output)
}




## Pablo: use this?
## Named matrix of thetas -> all of the model and predicted probs
## Recall these thetas have theta_i,j: effect of j on i.
MHN_from_thetas <- function(thetas) {
    oindex <- order(colnames(thetas))
    thetas <- thetas[oindex, oindex]
    output <- list()
    output[["MHN_theta"]] <- thetas
    output[["MHN_trans_rate_mat"]] <-
        theta_to_trans_rate_3_SM(thetas,
                                 inner_transition = inner_transitionRate_3_1)
    output[["MHN_trans_mat"]] <-
        trans_rate_to_trans_mat(output[["MHN_trans_rate_mat"]],
                                method = "competingExponentials",
                                paranoidCheck = TRUE)
    output[["MHN_td_trans_mat"]] <-
        trans_rate_to_trans_mat(output[["MHN_trans_rate_mat"]],
                                method = "uniformization",
                                paranoidCheck = TRUE)
    output[["MHN_exp_theta"]] <- exp(thetas)
    # output[["MHN_exp_theta"]] <- exp(thetas)
    return(output)
}




## Pablo: use this?
## poset as adjacency matrix and vector of lambdas
## both named -> all of the cbn output
CBN_from_poset_lambdas <- function(poset, lambdas) {
    stopifnot(identical(sort(colnames(poset)), sort(names(lambdas))))
    poset_as_data_frame <- poset_2_data_frame(poset)
    output <- list()
    output[["CBN_model"]] <- CBN_model_from_edges_lambdas(poset_as_data_frame,
                                                       lambdas)
    output <- c(output,
                CBN_model_2_output(output[["CBN_model"]]))
    
    return(output)
}


## Pablo: or use this if you use a data frame
## data frame with "From", "To", "Edges" and lambdas -> all of the cbn output
CBN_model_2_output <- function(model) {
    ## extra level of nesting
    tmpo <- cpm2tm(list(edges = model))
    output <- list()
    output[["CBN_trans_rate_mat"]] <- tmpo[["weighted_fgraph"]]
    output[["CBN_trans_mat"]] <- tmpo[["trans_mat_genots"]]
    output[["CBN_td_trans_mat"]] <-
        trans_rate_to_trans_mat(tmpo[["weighted_fgraph"]],
                                method = "uniformization")
    return(output)
}

## Attach the lambdas to the edges data frame
CBN_model_from_edges_lambdas <- function(edges, lambdas) {
    edges[["rerun_lambda"]] <- vapply(edges[, "To"],
                                      function(x) lambdas[x], 0.0)
    return(edges)
}


## Convert a poset matrix to a data frame with the "edges" structure
##   Adds the "Root" component, if not present
poset_2_data_frame <- function(poset) {
    stopifnot(identical(colnames(poset), rownames(poset)))
    if (!("Root" %in% colnames(poset))) {
        ## Add WT and WT connections
        not_connected <- which(colSums(poset) == 0)
        adjm2 <- cbind(rep(0, nrow(poset)), poset)
        adjm2 <- rbind(rep(0, ncol(adjm2)),
                       adjm2)
        adjm2[1, not_connected + 1] <- 1
        colnames(adjm2) <- rownames(adjm2) <- c("Root", colnames(poset))
        poset <- adjm2
    }
    elist <- igraph::get.edgelist(igraph::graph_from_adjacency_matrix(adjm2))

    edges <- data.frame(From = elist[, 1],
                        To = elist[, 2],
                        edge = paste0(elist[, 1], " -> ", elist[, 2]))
    return(edges)
}






## Pablo: use this one?
## poset as adjac. matrix, vector of lambdas, parent_set -> full HESBCN output
HESBCN_from_poset_lambdas_relation_probs <- function(poset, lambdas,
                                                     hesbcn_probs) {
    stopifnot(identical(sort(colnames(poset)), sort(names(lambdas))))
    poset_as_data_frame <- poset_2_data_frame(poset)

    ## Assign type of relationship randomly to nodes with >= 2 parents
    num_parents <- table(poset_as_data_frame[, "To"])
    singles <- names(which(num_parents == 1))
    parent_set <- vector(mode = "character", length = length(num_parents))
    names(parent_set) <- names(num_parents)
    parent_set[singles] <- "Single"
    lmultiple <- length(num_parents) - sum(num_parents == 1)
    if (lmultiple > 0) {
        ps_values <- sample(names(hesbcn_probs), size = lmultiple,
                            prob = hesbcn_probs, replace = TRUE)
        parent_set[which(num_parents > 1)] <- ps_values
    }
    
    stopifnot(identical(sort(names(parent_set)),
                        sort(names(lambdas))))
    
    output <- list()
    output[["HESBCN_parent_set"]] <- parent_set
    output[["HESBCN_model"]] <-
        HESBCN_model_from_edges_lambdas_parent_set(poset_as_data_frame,
                                                   lambdas,
                                                   parent_set)
    output <- c(output, HESBCN_model_2_output(output[["HESBCN_model"]],
                                              output[["HESBCN_parent_set"]]))
    return(output)
}





## Pablo: or use this if you use a data frame and parent frame and parent set
## data frame with "From", "To", "Edges" and lambdas -> all of the cbn output
HESBCN_model_2_output <- function(model, parent_set) {
    tmpo <- cpm2tm(list(edges = model, parent_set = parent_set))
    output <- list()
    output[["HESBCN_trans_rate_mat"]] <- tmpo[["weighted_fgraph"]]
    output[["HESBCN_trans_mat"]] <- tmpo[["trans_mat_genots"]]
    output[["HESBCN_td_trans_mat"]] <-
        trans_rate_to_trans_mat(tmpo[["weighted_fgraph"]],
                                method = "uniformization")
    ## output[["HESBCN_predicted_genotype_freqs"]] <- 
    ##     probs_from_trm(tmpo[["weighted_fgraph"]])
    return(output)
}

## Attach the lambdas to the edges data frame
## Same for Relation.
HESBCN_model_from_edges_lambdas_parent_set <- function(edges, lambdas,
                                                       parent_set) {
    edges[["Lambdas"]] <- vapply(edges[, "To"],
                                 function(x) lambdas[x], 0.0)
        
    edges[["Relation"]] <- vapply(edges[, "To"],
                                  function(x) parent_set[x], "")
    return(edges)
}

## A random poset where each gene depends on at most one parent
OT_random_poset <- function(ngenes, graph_density) {
    poset0 <- mccbn::random_poset(ngenes, graph_density = graph_density)
    ## Leave only one parent
    cosum <- colSums(poset0)
    if (all(cosum <= 1)) return(poset0)
    for (co in which(cosum > 1)) {
        ones <- which(poset0[, co] == 1)
        the_one <- sample(ones, size = 1)
        poset0[, co] <- 0L
        poset0[the_one, co] <- 1L
    }
    return(poset0)
}


## Pablo call this?
## poset as adjacency matrix, weights, epos -> full output, as from evam
##   weights: do not have Root
##   poset: one for OT, so no column with two or more parents
OT_from_poset_weights_epos <- function(poset, weights, epos) {
    stopifnot(identical(sort(colnames(poset)), sort(names(weights))))
    stopifnot(colSums(poset) <= 1)
    poset_as_data_frame <- poset_2_data_frame(poset)
    output <- list()
    output[["OT_model"]] <- OT_model_from_edges_lambdas(poset_as_data_frame,
                                                        weights)
    output <- c(output, OT_model_2_output(output[["OT_model"]],
                                          epos))
    return(output)
}


OT_model_from_edges_lambdas <- function(edges, weights) {
    edges[["OT_edgeWeight"]] <- vapply(edges[, "To"],
                                       function(x) weights[x], 0.0)
    return(edges)
}


## Pablo call this?
## OT model and epos -> full output, as from evam
OT_model_2_output <- function(model, epos) {
    ## We need to go back to the DAG representation
    ## Different from CBN: we obtain the probs. of genotypes
    ## using a call in Oncotree, that expects and oncotree.fit object.
    
    if (any(model$OT_edgeWeight > 1) || any(model$OT_edgeWeight < 0))
        stop("OncoBN's thetas must be between 0 and 1.")
    
    tmpo <- cpm2tm(list(edges = model))
    output <- list()
    output[["OT_f_graph"]] <- tmpo[["weighted_fgraph"]]
    output[["OT_trans_mat"]] <- tmpo[["trans_mat_genots"]]
    output[["OT_eps"]] <- c(epos = epos, eneg = 0)
    otm <- OT_model_2_predict_genots(model,
                                     epos)
    output[["OT_predicted_genotype_freqs"]] <- otm[["predicted_genotype_freqs"]]
    output[["OT_fit"]] <- otm[["fit"]]
    return(output)
}

OT_model_2_predict_genots <- function(model, epos) {
    ## Minimal checks, to stop if called incorrectly
    ## such as from shiny app
    if ("Relation" %in% colnames(model))
        stop("OT does not take a Relation column")
    if (any(table(model$To) > 1))
        stop("In OT models each gene has at most one parent.")
    
    ## Obtain the adjacency matrix and ensure adjacency matrix
    ## and weights have genes in same order
    adjm <- igraph::as_adjacency_matrix(
                       igraph::graph_from_data_frame(model[, c("From", "To")]))
    stopifnot(colnames(adjm)[1] == "Root")
    stopifnot(colnames(adjm) == rownames(adjm))
    ## Sort column names
    cnadjm_nor <- sort(setdiff(colnames(adjm), "Root"))
    adjm <- adjm[c("Root", cnadjm_nor), c("Root", cnadjm_nor)]
    weights <- model$OT_edgeWeight
    names(weights) <- model$To
    weights <- weights[cnadjm_nor]
    stopifnot(colnames(adjm)[-1] == names(weights))

    otfit <- oncotree_fit_from_adjm_weights_epos(adjm = as.matrix(adjm),
                                                 weights = weights,
                                                 epos = epos)
    ## We allow for errors from the model, not observational errors
    ##   epos >= 0, but eneg = 0
    ##  We set edge.weights to estimated. Observed ones are set to NA
    ##  As we use with.errors = TRUE, argument to edge.weights is not needed
    ##  but used to be explicit.
    preds <- distribution.oncotree(otfit,
                                   with.probs = TRUE,
                                   with.errors = TRUE,
                                   edge.weights = "estimated")
    preds <- dist_oncotree_output_2_named_genotypes(preds)
    stopifnot(isTRUE(all.equal(sum(preds), 1)))
    return(list(predicted_genotype_freqs = preds,
                fit = otfit))
}


## Adjacency matrix, weights, epos error -> list like that from oncotree.fit
## Adjacency matrix contains Root, weights do not.
oncotree_fit_from_adjm_weights_epos <- function(adjm, weights, epos) {
    otf <- list()
    otf$data <- NA
    otf$nmut <- ncol(adjm) 
    otf$parent <- oncotree_fit_parent_from_adjm_weights(adjm, weights)
    otf$eps <- c(epos = epos, eneg = 0)
    return(otf)
}

## Adjacency matrix and weights -> list like that from parent component of
## oncotree.fit
##      Adjacency matrix contains Root, weights do not.
oncotree_fit_parent_from_adjm_weights <- function(adjm, weights) {
    child <- rep("ERROR", times = ncol(adjm))
    parent <- rep("ERROR", times = ncol(adjm))
    parent.num <- rep(-99, times = ncol(adjm))
    for (p in seq_len(ncol(adjm))) {
        child[p] <- colnames(adjm)[p]
        tmp <- which(adjm[, p] == 1)
        if (length(tmp) == 0) {
            parent.num[p] <- 0
            parent[p] <- ""
        } else if (length(tmp) == 1) {
            parent.num[p] <- unname(tmp)
            parent[p] <- names(tmp)
        } else {
            stop("More than one parent")
        }
    }
    return(list(child = child,
                parent = parent,
                parent.num = parent.num,
                obs.weight = NA,
                est.weight = c(1, weights)
                ))
}



## o1 <- oncotree_fit_from_dag_weights_epos(ab, runif(5), 0.1)


## ## this is the right call
## o1p <- distribution.oncotree(o1, with.probs = TRUE, with.errors = TRUE, edge.weights = "estimated")

## ## And note these are identical
## o2 <- o1
## o2$parent$obs.weight <- runif(length(o1$parent$est.weight))
## o2p <- distribution.oncotree(o2, with.probs = TRUE, with.errors = TRUE, edge.weights = "estimated")
## stopifnot(all.equal(o1p$Prob, o2p$Prob))


## distribution.oncotree(o1, with.probs = TRUE, with.errors = FALSE, edge.weights = "estimated")


## distribution.oncotree(otf, with.probs = TRUE, with.errors = FALSE, edge.weights = "estimated")
## distribution.oncotree(otf, with.probs = TRUE, with.errors = TRUE, edge.weights = "estimated")


## Pablo calls this?
OncoBN_from_poset_thetas_epsilon_model <- function(poset,
                                                thetas,
                                                epsilon,
                                                model) {
    stopifnot(identical(sort(colnames(poset)), sort(names(thetas))))
    poset_as_data_frame <- poset_2_data_frame(poset)

    ## Assign type of relationship to nodes with >= 2 parents
    num_parents <- table(poset_as_data_frame[, "To"])
    singles <- names(which(num_parents == 1))
    parent_set <- vector(mode = "character",
                         length = length(num_parents))
    names(parent_set) <- names(num_parents)
    parent_set[singles] <- "Single"
    lmultiple <- length(num_parents) - sum(num_parents == 1)
    if (lmultiple > 0) {
        if (model == "CBN") {
            ps_value <- "AND"
        } else if (model == "DBN") {
            ps_value <- "OR"
        } else {
            stop("No valid model")
        }
        parent_set[which(num_parents > 1)] <- ps_value
    }
    
    stopifnot(identical(sort(names(parent_set)),
                        sort(names(thetas))))

    output <- list()
    output[["OncoBN_parent_set"]] <- parent_set
    output[["OncoBN_fitted_model"]] <- model
    output[["OncoBN_likelihood"]] <- NA
    output[["OncoBN_model"]] <-
        OncoBN_model_from_edges_thetas_parent_set(poset_as_data_frame,
                                                   thetas,
                                                   parent_set)

    output <- c(output, OncoBN_model_2_output(output[["OncoBN_model"]],
                                              epsilon))
    return(output)
}



OncoBN_model_from_edges_thetas_parent_set <- function(edges,
                                                       thetas, parent_set) {
    edges[["theta"]] <- vapply(edges[, "To"],
                                       function(x) thetas[x], 0.0)
    edges[["Relation"]] <- vapply(edges[, "To"],
                                  function(x) parent_set[x], "")
    return(edges)
}


## Pablo
OncoBN_model_2_output <- function(model, epsilon) {
    ## We need to go back to the DAG representation
 
    if (any(model$theta > 1) || any(model$theta < 0))
        stop("OncoBN's thetas must be between 0 and 1.")
    
    tmpo <- cpm2tm(list(edges = model))
    output <- list()
    output[["OncoBN_f_graph"]] <- tmpo[["weighted_fgraph"]]
    output[["OncoBN_trans_mat"]] <- tmpo[["trans_mat_genots"]]
    output[["OncoBN_epsilon"]] <- epsilon
    output[["OncoBN_predicted_genotype_freqs"]] <-
        OncoBN_model_2_predict_genots(model,
                                      epsilon = epsilon)
    return(output)
}



OncoBN_model_2_predict_genots <- function(model, epsilon) {
    ## Create a representation as used by OncoBN
    ## We need components: graph, theta, model, epsilon
    ## edgelist and score not added.

    ## Minimal checks, to stop if called incorrectly
    ## such as from shiny app
    if (any(model$Relation == "XOR"))
        stop("OncoBN does not accept XOR relations.")
    if (sum(c("OR", "AND") %in% model$Relation ) > 1)
        stop("OncoBN does not accept, in the same model, ",
             "both AND and OR relations")
    
    obnfit <- list()
    ## In OncoBN what we call Root is called WT
    model$From[model$From == "Root"] <- "WT"

    thetadf <- aggregate(theta ~ To, data = model, FUN = unique)
    thetav <- thetadf$theta
    names(thetav) <- thetadf$To
    theta <- thetav[sort(names(thetav))]
    names_g <- names(theta)
    obnfit[["theta"]] <- theta
    
    ## Make sure we have a graph
    ## that is always consistently ordered to prevent
    ## https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1048814030
    am <- igraph::as_adjacency_matrix(
                      igraph::graph_from_data_frame(model[, c("From", "To")]))
    am <- am[c("WT", names_g), c("WT", names_g)]
    obnfit[["graph"]] <- igraph::graph_from_adjacency_matrix(am)
    
    obnfit[["model"]] <- ifelse(any(model$Relation == "AND"), "CBN", "DBN")
    obnfit[["epsilon"]] <- epsilon
    pred_genots <- DBN_prob_genotypes(obnfit, sort(names(thetav)))
    pred_genots <- DBN_est_genots_2_named_genotypes(pred_genots)

    return(pred_genots)
}


## Given a data frame with columns
## exactly equal to either c(Genotype, Counts) or c(Genotype, Freq)
## where Genotype is a string
## or a named vector
## return a data set subjects as rows and columns as genes,
## with 0/1
counts_to_data_no_e <- function(x) {
    if (is.vector(x)) {
        stopifnot(!is.null(names(x)))
    } else if (is.data.frame(x)) {
        stopifnot(isTRUE(
            identical(names(x), c("Genotype", "Counts"))
            ||
            identical(names(x), c("Genotype", "Freq"))
        ))
        ## As per previous check we know x[, 2] is Counts or Freq
        if (is.data.frame(x)) x <- stats::setNames(x[, 2],
                                                   x$Genotype)
    } else {
        stop("Input must be a data frame or a named vector")
    }

    genes <- setdiff(unique(unlist(strsplit(names(x), ", "))),
                     "WT")
    ## Just in case
    genotypes <- canonicalize_genotype_names(names(x))

    ngenes <- length(genes)
    genotbin <- lapply(strsplit(genotypes, ", "),
                       function(u) {
                           wg <- which(genes %in% u)
                           v <- rep(0, ngenes)
                           v[wg] <- 1
                           return(v)})
    genotbin <- do.call(rbind, genotbin)
    colnames(genotbin) <- genes
    rr <- rep(1:nrow(genotbin), x)
    return(genotbin[rr, , drop = FALSE])
}



## From the identically named function in OncoSimulR
add_noise <- function(x, properr) {
    stopifnot(is.matrix(x))
    if (properr <= 0) {
        return(x)
    }
    else {
        if (properr > 1)
            stop("Proportion with error cannot be > 1")
        nn <- prod(dim(x))
        flipped <- sample(nn, round(nn * properr))
        x[flipped] <- as.integer(!x[flipped])
        return(x)
    }
}


## Given named vector or data frame with columns exactly
## equal to either c(Genotype, Counts) or c(Genotype, Freq)
## return a data set subjects as rows and columns as genes,
## with 0/1.
## e: noise error, as fraction (i.e., 0 to 1)
genotypeCounts_to_data <- function(x, e) {
    d <- counts_to_data_no_e(x)
    if (e > 0) d <- add_noise(d, e)
    return(d)
}

## The revert of counts to data. Return as standard order
## Similar to OncoSimulR's sampledGenotypes
## but using genot_matrix_2_vector.
## data: 0/1 data, in a data.frame or matrix
## out: one of vector or data.frame
## omit_0: if genotypes with 0 counts should be omitted
data_to_counts <- function(data, out,
                           omit_0 = FALSE) {
    stopifnot(out %in% c("vector", "data.frame"))

    if (is.data.frame(data)) data <- data.matrix(data)

    stopifnot(!is.null(colnames(data)))
    stopifnot(length(colnames(data)) == ncol(data))
    stopifnot(isTRUE(all(data %in% c(0, 1))))
    t_genots_string <- table(genot_matrix_2_vector(data))
    v_genots_string <- as.vector(t_genots_string)
    names(v_genots_string) <- names(t_genots_string)

    o_genots_string <- reorder_to_standard_order(v_genots_string)

    if (omit_0) {
        o_genots_string <- na.omit(o_genots_string)
        attributes(o_genots_string)$na.action <- NULL
    } else {
        o_genots_string[is.na(o_genots_string)] <- 0
    }
    stopifnot(nrow(data) == sum(o_genots_string))

    
    if (out == "vector") {
            return(o_genots_string)
    } else if (out == "data.frame") {
        df <- data.frame(Genotype = names(o_genots_string),
                         Counts   = o_genots_string)
        rownames(df) <- seq_len(nrow(df))
        return(df)
    } else {
        stop("Incorrect out option")
    }
}



## x: named vector of genotype probabilities
## ngenes: number of genes
## gene_names: gene names
## N: size of the sample
## output: ordered sample. Ordered as in pD and MHN

## A common construction that is used in many places.
genot_probs_2_pD_ordered_sample <- function(x, ngenes, gene_names, N, out) {
    this_sample <- sample(names(x), size = N, prob = x, replace = TRUE)
    return(
        sample_to_named_pD_ordered_out(the_sample = this_sample,
                                          ngenes = ngenes,
                                          gene_names = gene_names,
                                          out = out)
    )
}


## wrap the call to sample_to_pD_order so we get a named vector
sample_to_named_pD_ordered_out <- function(the_sample, ngenes, gene_names,
                                        out = c("vector", "data.frame")) {
    out <- match.arg(out)
    counts <- sample_to_pD_order(the_sample, ngenes, gene_names)
    pD_sorted_genotypes <- generate_pD_sorted_genotypes(ngenes,
                                                        gene_names)
    
    ordered_vector <- stats::setNames(counts, pD_sorted_genotypes)
    
    if (out == "vector") {
        return(ordered_vector)
    } else if (out == "data.frame") {
        return(data.frame(Genotype = names(ordered_vector),
                          Counts = ordered_vector))
    } else {
        stop("Incorrect out option")
    }
}



