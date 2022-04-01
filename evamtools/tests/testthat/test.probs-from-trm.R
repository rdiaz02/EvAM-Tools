## The testing is convoluted, but the original code
## in MHN does not compute the probabilities from a transition rate
## matrix but from the thetas which, in addition, can
## be an object without row/column names.

t1 <- Sys.time()

## Given some gene names, generate some random theta and from it
## the vector of probabilities. Return the theta and the probs.
## Since it is easy to make mistakes on what numbers correspond
##  to what genotype we double check sorting code.
probs_from_Schill <- function(gene_names) {
    ## Make sure we return named vector
    ## But do not order thetas by gene names, so as
    ## to also catch possible errors with sorting
    n <- length(gene_names)
    thetas <- Random.Theta(n = n)
    rownames(thetas) <- colnames(thetas) <- gene_names
    timep <- system.time(p <- Generate.pTh(thetas))["elapsed"]
    names(p) <- generate_pD_sorted_genotypes(n, gene_names,
                                                      sort_gene_names = FALSE)

    ## What if we had sorted?
    thetas2 <- thetas
    rownames(thetas2) <- colnames(thetas2) <- gene_names
    oindex <- evam_string_order(colnames(thetas2))
    thetas2 <- thetas2[oindex, oindex]
    rm(oindex)
    p2 <- Generate.pTh(thetas2)
    names(p2) <- generate_pD_sorted_genotypes(n, gene_names,
                                                       sort_gene_names = TRUE)

    ## gn_p_can might be the same as names(p) if gene_names was ordered
    gn_p_can <- canonicalize_genotype_names(names(p))
    ## gn_p2 should be the same
    gn_p2_can <- canonicalize_genotype_names(names(p2))
    stopifnot(identical(gn_p2_can, names(p2)))
    rm(gn_p2_can)
    
    p_can <- p
    names(p_can) <- gn_p_can
    ## We have only renamed
    stopifnot(isTRUE(all(p_can == p)))
    
    p2_as_p <- p2[names(p_can)]
    stopifnot(all.equal(p2_as_p, p_can))

    ## message("Time for Generate.pTh = ", timep)
    return(list(
        thetas = thetas,
        p = p,
        p_genes_sorted_in_genotypes = p_can,
        p_genes_and_genotypes_sorted = p2))
}


probs_from_theta_evam <- function(theta) {
    oindex <- evam_string_order(colnames(theta))
    theta <- theta[oindex, oindex]

    ## The previous sorting ensures genes are sorted in genotypes;
    ## genotypes themselves not necessarily sorted
    ## as in pD
    trm <- theta_to_trans_rate_3_SM(theta,
                                    inner_transition = inner_transitionRate_3_1)
    tptm <- system.time(p <- probs_from_trm(trm))["elapsed"]
    ## message("time in probs_from_trm = ", tptm)
    return(list(
        p = p,
        trm = trm
    ))
}


compare_schill_evam_probs <- function(gene_names) {
    s1 <- probs_from_Schill(gene_names)
    pe <- probs_from_theta_evam(s1$thetas)
    pe_as_schill <- pe$p[names(s1$p_genes_sorted_in_genotypes)]
    expect_equal(pe_as_schill, s1$p_genes_sorted_in_genotypes)
}


weird_gene_names <- function(ngenes) {
    gn <- vector(mode = "character", length = ngenes)
    for(g in seq_len(ngenes)) {
        l1 <- sample(c(LETTERS, letters), 1)
        rest <-  c(sample(
            c(
                sample(c("_",  "=", "?", "#", "@", "%", "&", "!"), 4, replace = TRUE)
              , sample(0:9, 4, replace = TRUE)
              , sample(c(letters, LETTERS), 4, replace = TRUE)
            )))
        
        rest <- paste(rest, sep = "", collapse = "")
        gn[g] <- paste0(l1, rest)
    }
    return(gn)
}



test_that("Same results as original MHN code by Schill", {
    for(i in 1:10) {
        ng <- 5
        compare_schill_evam_probs(LETTERS[1:ng])
        compare_schill_evam_probs(sample(LETTERS[1:ng]))
    }

    
    ## ## Overkill and slow. Do a single one with 10 genes
    compare_schill_evam_probs(sample(LETTERS[1:10]))
    ## for(i in 1:10) {
    ##     ng <- 11
    ##     compare_schill_evam_probs(LETTERS[1:ng])
    ##     compare_schill_evam_probs(sample(LETTERS[1:ng]))
    ## }
    
    for(i in 1:5) {
        ng <- sample(2:6, 1)
        wg <- weird_gene_names(ng)
        compare_schill_evam_probs(wg)
        compare_schill_evam_probs(wg)
    }
})


test_that("population_sample_from_trm: pre and not precompute", {
    for (i in 1:2) {
        ng <- sample(3:4, size = 1)
        thetas <- Random.Theta(n = ng)
        rownames(thetas) <- colnames(thetas) <- LETTERS[1:ng]
        trm <- theta_to_trans_rate_3_SM(thetas,
                                        inner_transition = inner_transitionRate_3_1)
        seed <- round(runif(1, 1, 1e8))
        set.seed(seed)
        ## Crucial not to use multiple cores or random numbers differ
        p1 <- population_sample_from_trm(trm, 1e3, rep(1, 1e3), pre_compute = TRUE,
                                         cores = 1)
        set.seed(seed)
        p2 <- population_sample_from_trm(trm, 1e3, rep(1, 1e3), pre_compute = FALSE,
                                         cores = 1)
        expect_equal(p1, p2)
        ## cat("\n seed = ", seed, "\n")
        set.seed(NULL)
        }
})



if(FALSE) {
    test_that("Compare the probs obtained with those from sampling", {
        ## Going down this route is not worth it. This is just an example.
        ## But these tests would take a long time to code and run for them
        ## to test things well automatically.
        ## See the tests for comparing between the sampling and the
        ## sampling from MC-CBN and MHN

        data(every_which_way_data)
        Dat1 <- every_which_way_data[[16]][1:40, 2:6]
        out <- suppressMessages(evam(Dat1,
                                     methods = c("CBN", "MHN", "HESBCN")))

        p_cpm <- probs_from_trm(out$CBN_trans_rate_mat, all_genotypes = TRUE)
        p_cpm <- reorder_to_pD(p_cpm)
        
        ## slows down a lot with larger n
        p_cpm_sampl <- population_sample_from_trm(out$CBN_trans_rate_mat,
                                                              n_samples = 1e7)
        
        v_sampl <- sample_to_pD_order(p_cpm_sampl$obs_events,
                                                  ngenes = ncol(out$analyzed_data),
                                                  gene_names = colnames(out$analyzed_data))
        
        names(v_sampl) <- generate_pD_sorted_genotypes(n_genes = ncol(out$analyzed_data),
                                                                gene_names = colnames(out$analyzed_data))
        v_sampl <- v_sampl/sum(v_sampl)

        t(apply(cbind(p_cpm, v_sampl), 1,
                function(y) {return(c(diff = y[1] - y[2],
                                      rel_diff = abs(y[1]- y[2])/min(y)))}))

        ## ######################################################################
        ## HESBCN. Examples where genotype with all loci mutated is not accessible


        set.seed(2)
        d1 <- data.frame(A = sample(c(1, 0), prob = c(0.7, 0.2), size = 200, replace = TRUE),
                         B = sample(c(1, 0), prob = c(0.85, 0.2), size = 200, replace = TRUE))
        d1$C <- 0
        d1$D <- 0
        d1$E <- 0
        d1$F <- 0
        d1$C[(d1$A == 1) & (d1$B == 1)] <- 1
        d1$D[(d1$A == 1) | (d1$B == 1)] <- 1
        d1$D[100:200] <- 0
        d1$E[xor((d1$A == 1), (d1$B == 1))] <- 1
        d1$F[(d1$C == 1) & (d1$D == 1)] <- 1
        d2 <- rbind(d1,
                    data.frame(A = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                               B = sample(c(1, 0), size = 25, prob = c(0.5, 0.2), replace = TRUE),
                               C = 0,
                               D = sample(c(1, 0), size = 25, replace = TRUE),
                               E = 0,
                               F = 0))
        d3 <- d2
        d3$C[(d3$A == 1) & (d3$B == 1)] <- 1

        ## Examples that mix output
        ## d3_2 <- do_HESBCN(d3, seed = 31)  ## AND, XOR, Single
        set.seed(NULL)

        out2 <- evam(d3, methods = c("HESBCN"), hesbcn_opts = list(seed = 31))
        
        p2_cpm <- probs_from_trm(out2$HESBCN_trans_rate_mat, all_genotypes = TRUE)
        p2_cpm <- reorder_to_pD(p2_cpm)
        
        ## slows down a lot with larger n
        p2_cpm_sampl <- population_sample_from_trm(out2$HESBCN_trans_rate_mat,
                                                               n_samples = 1e6)
        
        v2_sampl <- sample_to_pD_order(p2_cpm_sampl$obs_events,
                                                   ngenes = ncol(out2$analyzed_data),
                                                   gene_names = colnames(out2$analyzed_data))
        
        names(v2_sampl) <- generate_pD_sorted_genotypes(n_genes = ncol(out2$analyzed_data),
                                                                 gene_names = colnames(out2$analyzed_data))

        chisq.test(x = v2_sampl[p2_cpm > 0], p = p2_cpm[p2_cpm > 0], B = 2000)
        
        v2_sampl <- v2_sampl/sum(v2_sampl)

        comps <- t(apply(cbind(p2_cpm, v2_sampl), 1,
                         function(y) {return(c(diff = y[1] - y[2],
                                               rel_diff = abs(y[1]- y[2])/min(y)))}))
        ## So largest relative difference is 3% (of the smallest of the two)
        ## And largest abs. diff. is 6e-4 so 600 observations
        summary(comps0)

        ## Run another, much smaller
        p2b_cpm_sampl <- population_sample_from_trm(out2$HESBCN_trans_rate_mat,
                                                                n_samples = 1e5)
        
        v2b_sampl <- sample_to_pD_order(p2b_cpm_sampl$obs_events,
                                                    ngenes = ncol(out2$analyzed_data),
                                                    gene_names = colnames(out2$analyzed_data))
        
        names(v2b_sampl) <- generate_pD_sorted_genotypes(n_genes = ncol(out2$analyzed_data),
                                                                  gene_names = colnames(out2$analyzed_data))


        chisq.test(x = v2b_sampl[p2_cpm > 0], p = p2_cpm[p2_cpm > 0], B = 2000)
        
        
        v2b_sampl <- v2b_sampl/sum(v2b_sampl)

        comps <- t(apply(cbind(p2_cpm, v2b_sampl), 1,
                         function(y) {return(c(diff = y[1] - y[2],
                                               rel_diff = abs(y[1]- y[2])/min(y)))}))
        summary(comps)
    })
}

cat("\n Done test.probs-from-trm.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
