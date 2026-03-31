t1 <- Sys.time()


test_that("Genotype frequency predictions same as from OncoTree" , {
    require(OncoBN)
    require(Oncotree)
    data(ov.cgh)
    data(examples_csd)
    set.seed(NULL)

    ## Create our simple data frame with From, To, and probs/weights
    ## from the fitted object. We ARE NOT copying the predicted
    ## genotype freqs, just the model
    ## This is coming from do_OncoBN
    build_manual_obnfit <- function(fit, gene_names) {
        model <- fit$model
        thetas <- fit$theta
        adjm <- igraph::as_adjacency_matrix(
                            igraph::make_directed_graph(fit$edgelist))
        gn <- gene_names
        adjm <- adjm[c("WT", gn), c("WT", gn)]
        new_graph <- igraph::graph_from_adjacency_matrix(adjm)
        fit$graph <- new_graph

        dbn_out <- igraph::as_data_frame(fit$graph)
        colnames(dbn_out) <- c("From", "To")
        dbn_out$From[dbn_out$From == "WT"] <- "Root"
        dbn_out$edge <- mapply(function(x, y) { paste0(x, " -> ", y) },
                               dbn_out$From, dbn_out$To)
        dbn_out$theta <- thetas[dbn_out$To]

        dbn_out$Relation <- ifelse(model == "DBN", "OR", "AND")
        dbn_out$Relation[dbn_out$From == "Root"] <- "Single"
        return(dbn_out)
    }


    ## I modify this evamtools function, which is just
    ## a wrapper to GA_Likelihood. I need the wrapper
    ## because of the still open issue
    ## https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1033260046
    ## I modify it, w.r.t. evamtools one, by not passing gene_names:
    ## it gets them from the fitted
    ## model. But we must ensure the thetas are ordered
    ## correctly too.
    ## This, by the way, is another check we are doing the
    ## fixing of the problems of the G matrix and the thetas order
    ## properly (as we do it differently from what is done in evamtools'
    ## DBN_prob_genotypes
    DBN_prob_genotypes_mod <- function(fit, epsilon = 0) {

        ## No longer needed, as Lik.genotype does the right thing now.
        ## Not quite: it doesn't and I am moving out common ops.
        ## Pre-check same order and proper naming
        ## See https://github.com/phillipnicol/OncoBN/issues/3#issuecomment-1048814030
        G <- t(as.matrix(igraph::as_adjacency_matrix(fit$graph)))

        n <- ncol(G) - 1
        genotypes <- expand.grid(replicate(n, 0:1, simplify = FALSE))
        colnames(genotypes) <- colnames(G)[-1]

        stopifnot(colnames(G)[1] == "WT")

        thetas <- fit$theta[colnames(G)[-1]]
        stopifnot(colnames(G)[-1] == names(thetas))

        ## ## Extract that code, and call GA_Likelihood.
        ## genotypes$Prob <- apply(genotypes, 1,
        ##                         function(x) old_Lik.genotype(fit, x))

        ## Faster
        genotypes2 <- cbind(1, genotypes)
        ## theta.in <- c(1, fit$theta)
        theta.in <- c(1, thetas)
        ## epsilon the third argument to the function; o.w. default is fit$epsilon
        genotypes$Prob <- apply(genotypes2, 1,
                                function(x) OncoBN:::GA_Likelihood(x, G,
                                                                   theta.in,
                                                                   epsilon,
                                                                   fit$model))
        return(genotypes)
    }

    ## From the prediction matrix, return a named vector
    named_pred_vector <- function(x) {
        gene_names <- colnames(x[, -ncol(x)])
        names_genots <- apply(x[, -ncol(x)], 1,
                              function(v) paste(evamtools:::evam_string_sort(gene_names[v == 1])))
        names_genots <- unlist(lapply(names_genots, function(u) paste(u, collapse = ", ")))
        names_genots <- unname(names_genots)
        names_genots[names_genots == ""] <- "WT"
        preds <- x$Prob
        names(preds) <- names_genots
        return(preds)
    }

    ## Doing one should be enough, but we do a few rounds.
    ## The first uses the very original data set in Oncotree
    ## The rest, change the gene names and some data radomly
    iter <- 12
    for (i in 0:iter) {
        if (i == 0) df1 <- ov.cgh
        if (i == 1) df1 <- examples_csd$csd$AND$data
        if (i == 2) df1 <- examples_csd$csd$Linear$data
        if (i == 3) df1 <- examples_csd$csd$OR$data
        if (i == 4) df1 <- examples_csd$csd$XOR$data
        if (i == 5) df1 <- examples_csd$csd$c1$data
        if (i == 6) df1 <- examples_csd$csd$c4c2$data
        if (i > 6) {
            df1 <- ov.cgh
            colnames(df1) <- sample(LETTERS, ncol(df1))
            which_flip <- rbinom(prod(dim(df1)), 1, 0.05)
            which_flip <- matrix(which_flip, nrow = nrow(df1))
            df1 <- abs(df1 - which_flip)
        }


        obn_c <- fitCPN(df1, model = "CBN")
        obn_d <- fitCPN(df1, model = "DBN")

        obn_c_preds <- DBN_prob_genotypes_mod(obn_c, 0)
        obn_d_preds <- DBN_prob_genotypes_mod(obn_d, 0)

        obn_c_manual <- build_manual_obnfit(obn_c, colnames(df1))
        obn_d_manual <- build_manual_obnfit(obn_d, colnames(df1))

        obn_c_manual_preds <-
            evamtools:::OncoBN_model_2_output(obn_c_manual,
                                              epsilon = 0)$OncoBN_predicted_genotype_freqs
        obn_d_manual_preds <-
            evamtools:::OncoBN_model_2_output(obn_d_manual,
                                              epsilon = 0)$OncoBN_predicted_genotype_freqs

        obn_c_manual_preds <- obn_c_manual_preds[obn_c_manual_preds > 0]
        obn_d_manual_preds <- obn_d_manual_preds[obn_d_manual_preds > 0]

        obn_vector_preds_c <- named_pred_vector(obn_c_preds)
        obn_vector_preds_d <- named_pred_vector(obn_d_preds)

        obn_vector_preds_c <- obn_vector_preds_c[obn_vector_preds_c > 0]
        obn_vector_preds_d <- obn_vector_preds_d[obn_vector_preds_d > 0]

        ## Same genotypes
        expect_true(isTRUE(all.equal(sort(names(obn_c_manual_preds)),
                                     sort(names(obn_vector_preds_c)))))
        expect_true(isTRUE(all.equal(sort(names(obn_d_manual_preds)),
                                     sort(names(obn_vector_preds_d)))))

        obn_vector_preds_c <- obn_vector_preds_c[order(obn_vector_preds_c)]
        obn_c_manual_preds <- obn_c_manual_preds[names(obn_vector_preds_c)]

        obn_vector_preds_d <- obn_vector_preds_d[order(obn_vector_preds_d)]
        obn_d_manual_preds <- obn_d_manual_preds[names(obn_vector_preds_d)]

        expect_true(isTRUE(all.equal(obn_vector_preds_c, obn_c_manual_preds)))
        expect_true(isTRUE(all.equal(obn_vector_preds_d, obn_d_manual_preds)))


        ## More consistency checks.
        ## Use the usual OncoBN, allowing for epsilon
        ev_c <- evamtools:::do_OncoBN(df1, model = "CBN")
        ev_d <- evamtools:::do_OncoBN(df1, model = "DBN")
        ev_c_p <- ev_c$predicted_genotype_freqs
        ev_d_p <- ev_d$predicted_genotype_freqs

        ## The fitted by OncoBN and the hand DBN_prob_genotypes is
        ## the same if we use the epsilon from the fit
        obn_c_preds_epsilon <- named_pred_vector(DBN_prob_genotypes_mod(obn_c, obn_c$epsilon))
        obn_d_preds_epsilon <- named_pred_vector(DBN_prob_genotypes_mod(obn_d, obn_d$epsilon))

        expect_true(isTRUE(all.equal(ev_c_p, obn_c_preds_epsilon[names(ev_c_p)])))
        expect_true(isTRUE(all.equal(ev_d_p, obn_d_preds_epsilon[names(ev_d_p)])))

        ## When comparing the manual preds. with 0 epsilon
        ## against evamtools fits, things differ

        expect_false(isTRUE(all.equal(ev_c_p[names(obn_c_manual_preds)], obn_c_manual_preds)))
        expect_false(isTRUE(all.equal(ev_d_p[names(obn_d_manual_preds)], obn_d_manual_preds)))
        expect_false(isTRUE(all.equal(ev_c_p, obn_c_manual_preds[names(ev_c_p)])))
        expect_false(isTRUE(all.equal(ev_d_p, obn_d_manual_preds[names(ev_d_p)])))

        ## But the manual with the non-zero epsilon is identical to
        ## the evamtools with epsilon from fit, and thus (as per comparison
        ## above) to the predictions from OncoBN
        obn_c_manual_preds_eps <-
            evamtools:::OncoBN_model_2_output(obn_c_manual,
                                              epsilon = obn_c$epsilon)$OncoBN_predicted_genotype_freqs
        obn_d_manual_preds_eps <-
            evamtools:::OncoBN_model_2_output(obn_d_manual,
                                              epsilon = obn_d$epsilon)$OncoBN_predicted_genotype_freqs

        expect_true(isTRUE(all.equal(ev_c_p, obn_c_manual_preds_eps[names(ev_c_p)])))
        expect_true(isTRUE(all.equal(ev_d_p, obn_d_manual_preds_eps[names(ev_d_p)])))

    }
}
)

set.seed(NULL)

cat("\n Done test.OncoBN-predicted-genotype-freqs.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
