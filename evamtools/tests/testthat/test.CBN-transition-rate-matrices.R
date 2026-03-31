t1 <- Sys.time()

## Direct tests of CBN transition rate matrix entries against
## hand-computed values.
##
## For CBN (conjunctive Bayesian network / poset-based model), a gene
## can mutate only when all its parents are already mutated, at a rate
## equal to its lambda.  The weighted_fgraph entry from genotype g to
## genotype g' is the lambda of the single gene added (g' = g + {x}),
## provided all parents of x are in g; otherwise the entry is 0.
##
## This test checks:
##   (1) The set of accessible genotypes.
##   (2) Every non-zero entry by name, against the hand-computed value.
##   (3) The total count and sum of non-zero entries.
##
## Examples ex0-ex5 use the same edge sets as test.trans-rates-f-graphs.R.
## There they are tested indirectly (via path probabilities and
## comparison with cpm_access_genots_paths_w). Here we verify the
## actual rate-matrix entries directly.

test_that("CBN: weighted_fgraph entries against hand-computed values", {

    compare_cbn_trm <- function(ex, ref_genots, ref_entries) {
        wfg <- as.matrix(cpm2tm(ex)$weighted_fgraph)
        ## (1) Accessible genotypes
        expect_equal(sort(rownames(wfg)), sort(ref_genots))
        ## (2) Every expected entry
        for (nm in names(ref_entries)) {
            parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
            expect_equal(
                wfg[parts[1], parts[2]], ref_entries[[nm]],
                label = paste0("wfg[", parts[1], " | ", parts[2], "]")
            )
        }
        ## (3) Total count and sum of non-zero entries
        expect_equal(sum(wfg != 0), length(ref_entries))
        expect_equal(sum(wfg),      sum(ref_entries))
    }


    ## ------------------------------------------------------------------
    ## ex0: linear chain  Root -> B(1) -> C(2) -> D(3) -> A(4)
    ## Accessible genotypes: WT, B, B+C, B+C+D, A+B+C+D
    ## ------------------------------------------------------------------
    ex0 <- list(edges = data.frame(
        From = c("Root", "B", "C", "D"),
        To   = c("B",    "C", "D", "A"),
        rerun_lambda = c(1, 2, 3, 4)
    ))
    compare_cbn_trm(
        ex0,
        ref_genots  = c("WT", "B", "B, C", "B, C, D", "A, B, C, D"),
        ref_entries = c(
            "WT|B"               = 1,
            "B|B, C"             = 2,
            "B, C|B, C, D"       = 3,
            "B, C, D|A, B, C, D" = 4
        )
    )


    ## ------------------------------------------------------------------
    ## ex1: Root->A(1), Root->B(2), [A AND B]->C(3), C->D(4)
    ## C requires both A and B (conjunction); D requires C.
    ## Accessible: WT, A, B, A+B, A+B+C, A+B+C+D
    ## ------------------------------------------------------------------
    ex1 <- list(edges = data.frame(
        From = c("Root", "Root", "A", "B", "C"),
        To   = c("A",    "B",    "C", "C", "D"),
        rerun_lambda = c(1, 2, 3, 3, 4)
    ))
    compare_cbn_trm(
        ex1,
        ref_genots  = c("WT", "A", "B", "A, B",
                        "A, B, C", "A, B, C, D"),
        ref_entries = c(
            "WT|A"               = 1,
            "WT|B"               = 2,
            "A|A, B"             = 2,
            "B|A, B"             = 1,
            "A, B|A, B, C"       = 3,
            "A, B, C|A, B, C, D" = 4
        )
    )


    ## ------------------------------------------------------------------
    ## ex2: Root->A(10), Root->B(11), A->C(12), B->D(14)
    ## Two independent branches (A->C and B->D); no conjunction.
    ## Accessible: WT, A, B, A+B, A+C, B+D, A+B+C, A+B+D, A+B+C+D
    ## ------------------------------------------------------------------
    ex2 <- list(edges = data.frame(
        From = c("Root", "Root", "A", "B"),
        To   = c("A",    "B",    "C", "D"),
        rerun_lambda = c(10, 11, 12, 14)
    ))
    compare_cbn_trm(
        ex2,
        ref_genots  = c("WT", "A", "B", "A, B", "A, C", "B, D",
                        "A, B, C", "A, B, D", "A, B, C, D"),
        ref_entries = c(
            "WT|A"               = 10,
            "WT|B"               = 11,
            "A|A, B"             = 11,
            "A|A, C"             = 12,
            "B|A, B"             = 10,
            "B|B, D"             = 14,
            "A, B|A, B, C"       = 12,
            "A, B|A, B, D"       = 14,
            "A, C|A, B, C"       = 11,
            "B, D|A, B, D"       = 10,
            "A, B, C|A, B, C, D" = 14,
            "A, B, D|A, B, C, D" = 12
        )
    )


    ## ------------------------------------------------------------------
    ## ex3: Root->A(2), Root->B(3), A->C(4), C->D(5)
    ## B is fully independent; C requires A; D requires C.
    ## Accessible: WT, A, B, A+B, A+C, A+B+C, A+C+D, A+B+C+D
    ## ------------------------------------------------------------------
    ex3 <- list(edges = data.frame(
        From = c("Root", "Root", "A", "C"),
        To   = c("A",    "B",    "C", "D"),
        rerun_lambda = c(2, 3, 4, 5)
    ))
    compare_cbn_trm(
        ex3,
        ref_genots  = c("WT", "A", "B", "A, B", "A, C",
                        "A, B, C", "A, C, D", "A, B, C, D"),
        ref_entries = c(
            "WT|A"               = 2,
            "WT|B"               = 3,
            "A|A, B"             = 3,
            "A|A, C"             = 4,
            "B|A, B"             = 2,
            "A, B|A, B, C"       = 4,
            "A, C|A, B, C"       = 3,
            "A, C|A, C, D"       = 5,
            "A, B, C|A, B, C, D" = 5,
            "A, C, D|A, B, C, D" = 3
        )
    )


    ## ------------------------------------------------------------------
    ## ex4: Root->A(1), Root->B(2), A->C(3), A->D(4)
    ## B is independent; C and D both require A (fan-out from A).
    ## Accessible: WT, A, B, A+B, A+C, A+D,
    ##             A+B+C, A+B+D, A+C+D, A+B+C+D
    ## ------------------------------------------------------------------
    ex4 <- list(edges = data.frame(
        From = c("Root", "Root", "A", "A"),
        To   = c("A",    "B",    "C", "D"),
        rerun_lambda = c(1, 2, 3, 4)
    ))
    compare_cbn_trm(
        ex4,
        ref_genots  = c("WT", "A", "B", "A, B", "A, C", "A, D",
                        "A, B, C", "A, B, D", "A, C, D", "A, B, C, D"),
        ref_entries = c(
            "WT|A"               = 1,
            "WT|B"               = 2,
            "A|A, B"             = 2,
            "A|A, C"             = 3,
            "A|A, D"             = 4,
            "B|A, B"             = 1,
            "A, B|A, B, C"       = 3,
            "A, B|A, B, D"       = 4,
            "A, C|A, B, C"       = 2,
            "A, C|A, C, D"       = 4,
            "A, D|A, B, D"       = 2,
            "A, D|A, C, D"       = 3,
            "A, B, C|A, B, C, D" = 4,
            "A, B, D|A, B, C, D" = 3,
            "A, C, D|A, B, C, D" = 2
        )
    )


    ## ------------------------------------------------------------------
    ## ex5: Root->A(1), Root->B(2), Root->C(3), [A AND B AND C]->D(4)
    ## Three-parent conjunction: D requires all of A, B, C.
    ## Accessible: WT, A, B, C, A+B, A+C, B+C, A+B+C, A+B+C+D
    ## ------------------------------------------------------------------
    ex5 <- list(edges = data.frame(
        From = c("Root", "Root", "Root", "A", "B", "C"),
        To   = c("A",    "B",    "C",    "D", "D", "D"),
        rerun_lambda = c(1, 2, 3, 4, 4, 4)
    ))
    compare_cbn_trm(
        ex5,
        ref_genots  = c("WT", "A", "B", "C",
                        "A, B", "A, C", "B, C",
                        "A, B, C", "A, B, C, D"),
        ref_entries = c(
            "WT|A"               = 1,
            "WT|B"               = 2,
            "WT|C"               = 3,
            "A|A, B"             = 2,
            "A|A, C"             = 3,
            "B|A, B"             = 1,
            "B|B, C"             = 3,
            "C|A, C"             = 1,
            "C|B, C"             = 2,
            "A, B|A, B, C"       = 3,
            "A, C|A, B, C"       = 2,
            "B, C|A, B, C"       = 1,
            "A, B, C|A, B, C, D" = 4
        )
    )

})


cat("\n Done test.CBN-transition-rate-matrices.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
