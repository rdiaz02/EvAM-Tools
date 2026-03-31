t1 <- Sys.time()

test_that("Genotype frequency predictions same as from OncoTree" , {
    require(Oncotree)
    data(ov.cgh)
    data(examples_csd)
    set.seed(NULL)

    ## Doing one should be enough, but we do a few rounds.
    ## The first uses the very original data set in Oncotree
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
            ## The rest, change the gene names and some data radomly
            df1 <- ov.cgh
            colnames(df1) <- sample(LETTERS, ncol(df1))
            which_flip <- rbinom(prod(dim(df1)), 1, 0.05)
            which_flip <- matrix(which_flip, nrow = nrow(df1))
            df1 <- abs(df1 - which_flip)
        }

        ov.tree <- oncotree.fit(df1)

        ## force epos and eneg to 0
        ov.tree.e0 <- ov.tree
        ov.tree.e0$eps <- c(epos = 0, eneg = 0)

        ## Predictions from Oncotree
        pred1 <- distribution.oncotree(ov.tree.e0,
                                       with.probs = TRUE,
                                       with.errors = TRUE,
                                       edge.weights = "estimated")

        pred1_no0 <- pred1[pred1$Prob > 0, ]

        ## Construct the model "by hand"
        manual.ov <-
            data.frame(From = ov.tree$parent$parent[-1],
                       To =   ov.tree$parent$child[-1],
                       OT_edgeWeight = ov.tree$parent$est.weight[-1])

        pred_manual.ov <- evamtools:::OT_model_2_output(manual.ov, epos = 0)$OT_predicted_genotype_freqs
        pred_manual.ov <- pred_manual.ov[pred_manual.ov > 0]


        gene_names <- colnames(pred1_no0[, -c(1, ncol(pred1_no0))])

        names_genots <- apply(pred1_no0[, -c(1, ncol(pred1_no0))], 1,
                              function(v) paste(evamtools:::evam_string_sort(gene_names[v == 1])))
        names_genots <- unlist(lapply(names_genots, function(u) paste(u, collapse = ", ")))
        names_genots <- unname(names_genots)
        names_genots[names_genots == ""] <- "WT"

        ot_preds <- pred1_no0$Prob
        names(ot_preds) <- names_genots

        ## Same genotypes
        expect_true(isTRUE(all.equal(sort(names(pred_manual.ov)),
                                     sort(names_genots))))

        ot_preds <- ot_preds[order(names_genots)]
        ev_preds <- pred_manual.ov[names(ot_preds)]

        ## cbind(ot_preds, ev_preds)

        expect_true(isTRUE(all.equal(ot_preds, ev_preds)))
    }
})

set.seed(NULL)

cat("\n Done test.OT-predicted-genotype-freqs.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
