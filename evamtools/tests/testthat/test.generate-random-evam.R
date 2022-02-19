test_that("Exercise random_evam with different options", {

    for (mm in c("MHN", "CBN", "HESBCN")) {
        ou <- paste0(mm, "_predicted_genotype_freqs")
        o1 <- generate_random_evam(ngenes = 5, model = mm)
        expect_true(!any(is.na(o1[[ou]])))
        }

    for (mm in c("MHN", "CBN", "HESBCN")) {
        ng <- sample(2:6, 1)
        gnam <- sample(LETTERS, size = ng)
        ou <- paste0(mm, "_predicted_genotype_freqs")
        o1 <- generate_random_evam(gene_names = gnam, model = mm,
                                   mhn_sparsity = 0.2,
                                   cbn_hesbcn_graph_density = 0.75,
                                   cbn_hesbcn_lambda_min = 1/4,
                                   cbn_hesbcn_lambda_max = 6,
                                   hesbcn_probs = c(AND = 0.1,
                                                    OR = 0.6,
                                                    XOR = 0.3)
                                   )
        expect_true(!any(is.na(o1[[ou]])))
    }

    
    ng <- sample(2:6, 1)
    gnam <- sample(LETTERS, size = ng)
    ou <- paste0("HESBCN", "_predicted_genotype_freqs")
    expect_error(o1 <- generate_random_evam(gene_names = gnam, model = "HESBCN",
                               mhn_sparsity = 0.2,
                               cbn_hesbcn_graph_density = 0.75,
                               cbn_hesbcn_lambda_min = 1/4,
                               cbn_hesbcn_lambda_max = 6,
                               hesbcn_probs = c(AND = 0.1,
                                                UR = 0.6,
                                                XOR = 0.3)
                               ),
                 'identical(sort(names(hesbcn_probs)), c("AND", "OR", "XOR")) is not TRUE',
                 fixed = TRUE)
    
})
