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
    
    while(TRUE) {
        ## qii <- diag[row] ## Or qs in Gotovos et al terminology
        ## Special case of genotype that does not transition to anything
        if(diag[row] == 0.0) break 
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
    #   , trans_table = NA ## I do not know what this is for
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
                            output = c("sampled_genotype_freqs",
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
    valid_output <- c("sampled_genotype_freqs",
                      "state_counts",
                      "obs_genotype_transitions")

    if (any(!(output %in% valid_output ))) stop("Incorrect output specified")
    if (length(output) == 0) stop("No output specified")
    
    ## Set up
    retval <- list()
    n_states <- 2^n_genes
    sorted_genotypes <- generate_sorted_genotypes(n_genes, gene_names)

    ## Calculate frequencies: genotype frequencies
    if ("sampled_genotype_freqs" %in% output) {
        counts_tmp <- sample_to_pD_order(sim$obs_events, n_genes, gene_names)
        frequencies <- data.frame(
            Genotype = sorted_genotypes,
            Counts = counts_tmp
        )
        rownames(frequencies) <- NULL
        retval$sampled_genotype_freqs <- frequencies
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
        if(!(all(apso[last - 1, ] == apso[last, ])))
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
        state_counts <- sample_to_pD_order(unlisted_trajectories,
                                           n_genes, gene_names)
        state_counts <- data.frame(
            Genotype = sorted_genotypes,
            Counts = state_counts
        )
        retval$state_counts <- state_counts
        ## FIXME: Paranoid check. Will remove later
        if("obs_genotype_transitions" %in% output) {
            cstt2 <- colSums(tt2)[-1]
            ## No WT
            statecounts_vector <- state_counts[-1, "Counts"]
            names(statecounts_vector) <- state_counts[-1, "Genotype"]
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


sample_CPMs <- function(cpm_output
                      , N
                      , methods = c("OT", "OncoBN",
                                    "CBN", "MCCBN",
                                    "MHN", "HESBCN")
                      , output = c("sampled_genotype_freqs")
                        ## "obs_genotype_transitions",
                        ## "state_counts")
                        ## , obs_genotype_transitions = TRUE
                        ) {
    ## And I have "Source" for a source data type for the web server
    retval <- list()

    output <- unique(output)
    valid_output <- c("sampled_genotype_freqs",
                      "obs_genotype_transitions",
                      "state_counts")
    not_valid_output <- which(!(output %in% valid_output))
    if (length(not_valid_output)) {
        warning("Output(s) ",
                paste(output[not_valid_output], sep = ", ", collapse = ", "),
                " not among the available output.",
                " Ignoring the invalid output.")
        output <- methods[-not_valid_output]
    }
    if (length(output) == 0) stop("No valid output given.")
    if (any(c("state_counts", "obs_genotype_transitions") %in% output)) {
        message("For the requested output we will need to simulate ",
                "from the transition rate matrix.")
    }
    
    gene_names <- colnames(cpm_output$analyzed_data)
    n_genes <- length(gene_names)

    for (method in methods) {
        if (method %in% c("OT", "OncoBN")) {
            if (method == "OT") {
                tmp_data <- cpm_output$OT_predicted_genotype_freqs
                genots <- tmp_data[2:(ncol(tmp_data) - 1)]
            } else if (method == "OncoBN") {
                tmp_data <- cpm_output$OncoBN_predicted_genotype_freqs
                genots <- tmp_data[1:(ncol(tmp_data) - 1)]
            }

            genots_2 <- unlist(apply(genots, 1,
                                     function(x) paste(names(genots)[x == 1],
                                                       collapse = ", ")))
            names(genots_2) <- NULL
            genots_2[genots_2 == ""] <- "WT"

            tmp_genotypes_sampled <- sample_to_pD_order(
                sample(genots_2, size = N,
                       prob = tmp_data$Prob, replace = TRUE),
                ngenes = n_genes, gene_names = gene_names)

            retval[[sprintf("%s_sampled_genotype_freqs", method)]] <-
                data.frame(
                    Genotype = generate_sorted_genotypes(n_genes, gene_names),
                    Counts = tmp_genotypes_sampled
                )
            ## The next one is NOT implicitly available.
            ##   see OT_transition_matrices.org
            retval[[sprintf("%s_obs_genotype_transitions", method)]] <- NULL
        } else {
            ## evam always returns the method_trans_rate_mat
            ## even if just with an NA
            if (any(c("state_counts", "obs_genotype_transitions") %in% output)) {
                ## Need to simulate from trm
                trm <- cpm_output[[sprintf("%s_trans_rate_mat", method)]]
                if ((length(trm) == 1) && is.na(trm)) {
                    retval[[sprintf("%s_obs_genotype_transitions", method)]] <- NULL

                    ogt <- cpm_output[[sprintf("%s_predicted_genotype_freqs", method)]]
                    if ((length(ogt) == 1) && is.na(ogt)) { 
                        retval[[sprintf("%s_sampled_genotype_freqs", method)]] <- NULL
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
                    retval[[sprintf("%s_sampled_genotype_freqs",
                                    method)]] <- NULL
                } else {
                    tmp_genotypes_sampled <-
                        sample_to_pD_order(sample(names(genots_pred),
                                                  size = N,
                                                  prob = genots_pred,
                                                  replace = TRUE),
                                           ngenes = n_genes,
                                           gene_names = gene_names)

                    retval[[sprintf("%s_sampled_genotype_freqs", method)]] <-
                        data.frame(
                            Genotype = generate_sorted_genotypes(n_genes, gene_names),
                            Counts = tmp_genotypes_sampled)
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

sample_to_pD_order <- function(x, ngenes, gene_names = NULL) {
    if(is.null(gene_names)) gene_names <- LETTERS[seq_len(ngenes)]
    stopifnot(ngenes == length(gene_names))
    x <- as.data.frame(table(x), stringsAsFactors = FALSE)
    gene_names <- sort(gene_names)
    genot_int <- x[, 1]
    genot_int <- gsub("^WT$", "", genot_int, fixed = "FALSE")

    genot_int <- vapply(genot_int,
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
##  i.e., in the same order that MHN and generate_sorted_genotypes
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
    sorted_genots <- generate_sorted_genotypes(length(gene_n),
                                               gene_names = gene_n)

    if(!all(genots_n %in% sorted_genots))
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

    if(!all(genots_n %in% sorted_genots))
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

## Yes, this is much slower, like two orders of magnitude,
## than Schill's Generate.pTh. Still, about 0.3 to 0.4 seconds
## for 11 genes, and most than 90% spent in the checks.
probs_from_trm <- function(x,
                           tolerance = 10 * sqrt(.Machine$double.eps),
                           all_genotypes = TRUE) {
    p0 <-  c(1, rep(0, nrow(x) - 1))

    if (Matrix::nnzero(tril(x)))
        stop("Lower triangular not 0. Is this transposed?")
    Q <- t(x)
    diag(Q) <- -1 * colSums(Q)

    ## Equation 4 in Schill et al. Thus
    ## (I - Q) * p = p0
    I_Q <- Matrix::Diagonal(nrow(Q)) - Q

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
        p <- p4
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
