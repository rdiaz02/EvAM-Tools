## The testing is convoluted, but the original code
## in MHN does not compute the probabilities from a transition rate
## matrix but from the thetas which, in addition, can
## be an object without row/column names. 


## Given some gene names, generate some random theta and from it
## the vector of probabilities. Return the theta and the probs.
## Since it is easy to make mistakes on what numbers correspond
##  to what genotype we double check sorting code.
probs_from_Schill <- function(gene_names) {
    ## Make sure we return named vector
    ## But do not order thetas by gene names, so as
    ## to also catch possible errors with sorting
    n <- length(gene_names)
    thetas <- evamtools:::Random.Theta(n = n)
    rownames(thetas) <- colnames(thetas) <- gene_names
    timep <- system.time(p <- evamtools:::Generate.pTh(thetas))["elapsed"]
    names(p) <- evamtools:::generate_sorted_genotypes(n, gene_names,
                                                      sort_gene_names = FALSE)

    ## What if we had sorted?
    thetas2 <- thetas
    rownames(thetas2) <- colnames(thetas2) <- gene_names
    oindex <- order(colnames(thetas2))
    thetas2 <- thetas2[oindex, oindex]
    rm(oindex)
    p2 <- evamtools:::Generate.pTh(thetas2)
    names(p2) <- evamtools:::generate_sorted_genotypes(n, gene_names,
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
    oindex <- order(colnames(theta))
    theta <- theta[oindex, oindex]

    ## The previous sorting ensures genes are sorted in genotypes;
    ## genotypes themselves not necessarily sorted
    ## as in pD
    trm <- evamtools:::theta_to_trans_rate_3_SM(theta,
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
        ng <- 8
        compare_schill_evam_probs(LETTERS[1:ng])
        compare_schill_evam_probs(sample(LETTERS[1:ng]))
    }

    ## Overkill and slow
    ## for(i in 1:10) {
    ##     ng <- 11
    ##     compare_schill_evam_probs(LETTERS[1:ng])
    ##     compare_schill_evam_probs(sample(LETTERS[1:ng]))
    ## }
       
    for(i in 1:10) {
        ng <- sample(2:7, 1)
        wg <- weird_gene_names(ng)
        compare_schill_evam_probs(wg)
        compare_schill_evam_probs(wg)
    }
})



data(every_which_way_data)
Dat1 <- every_which_way_data[[16]][1:40, 2:6]
out <- suppressMessages(evam(Dat1,
                             methods = c("CBN", "MHN", "HESBCN")))

